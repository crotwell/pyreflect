#!/usr/bin/env python3

import asyncio
import os
import pprint
import json
import math
import numpy
import obspy
from obspy.taup import TauPyModel
from obspy.core.trace import Stats
from obspy.core.utcdatetime import UTCDateTime

from pyreflect import earthmodel, distaz, momenttensor, specfile

# path to mgenkennett, relative to simple subdir (add extra ../ )
reflectivityPath = '../../../RandallReflectivity'
distDeg = 45
azimuth = 45
origin_depth = 0
reduceVel=8.0  # set below based on earliest arrival
offset=-30
# good idea to have P in phase list to set reduceVel value
phaseList = ["P", "S", "SSvmp", "SPvmp", "SSvms"]
runName = 'simple'
if not os.path.isdir(runName):
    os.mkdir(runName)

# calc travel times to estimate model depth and slowness values
taumodel = TauPyModel(model="ak135" )
radiusOfEarth = 6371 # for flat to spherical ray param conversion, should get from model

arrivals = taumodel.get_pierce_points(source_depth_in_km=origin_depth/1000,
                                  distance_in_degree=distDeg,
                                  phase_list=phaseList)
maxDepth = 0
minRayParam = arrivals[0].ray_param
maxRayParam = minRayParam
DEPTH_INDEX=3
earliestArrival = arrivals[0]
for a in arrivals:
    if earliestArrival.time > a.time:
        earliestArrival = a
    if minRayParam > a.ray_param:
        minRayParam = a.ray_param
    if maxRayParam < a.ray_param:
        maxRayParam = a.ray_param
    for p in a.pierce:
        if p[DEPTH_INDEX] > maxDepth:
            maxDepth = p[DEPTH_INDEX]
maxDepth = round(math.ceil(maxDepth + 10)) # little bit deeper
minRayParam = minRayParam/radiusOfEarth # need to be flat earth ray params
maxRayParam = maxRayParam/radiusOfEarth

# calc red velocity of earliest arrival (P wave usually)
reduceVel = distaz.DistAz.degreesToKilometers(distDeg) / earliestArrival.time

# load model file for PREM down to 85 km depth

model = earthmodel.EarthModel.loadAk135f(maxDepth)
model.name = runName
momtensorvals = {'m_dd': -0.136e+26, 'm_nn': -0.073e+26, 'm_ee': 0.209e+26, 'm_nd': 0.074e+26,
                           'm_ed': -1.262e+26, 'm_ne': -0.310e+26}
model.momentTensor = momtensorvals
model.distance = {
            "type": earthmodel.DIST_SINGLE,
            "distance": distaz.DistAz.degreesToKilometers(distDeg),
            "azimuth": azimuth
        }
model.frequency = {
            "min": 0.0,
            "max": 1.0,
            "nyquist":1.0,
            "numtimepoints": 4096
        }
model.slowness = {
    "lowcut": minRayParam/2,
    "lowpass": minRayParam*0.9,
    "highpass": maxRayParam*1.1,
    "highcut": maxRayParam*1.5,
    "controlfac": 1.0
    }

# pretty print all the model parameters as json,
# this should help if you need to modify things:
model.writeToJsonFile(os.path.join(runName, "model.json"))
model.writeToFile(os.path.join(runName, "simplemodel.ger"))
model.eft().writeToFile(os.path.join(runName, "simplemodel_eft.ger"))

async def mgenkennett(model):
    mgenkennett = f'{reflectivityPath}/mgenkennett'
    return await runCode(mgenkennett, model)


async def runCode(code, model):
    modelFilename = "model_mgen_eft.ger"
    model.eft().writeToFile(os.path.join(runName, modelFilename))
    proc = await asyncio.create_subprocess_shell(
        f"{code} {modelFilename}",
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=runName)

    stdout, stderr = await proc.communicate()

    print(f'[{code} exited with {proc.returncode}]')
    if stdout:
        print(f'[stdout]\n{stdout.decode()}')
    if stderr:
        print(f'[stderr]\n{stderr.decode()}')
    if proc.returncode != 0:
        print(f'Warning, {code} did not exit successfully')
    return proc.returncode

async def runSpec2zrt(reduceVel, offset):
    code = f'{reflectivityPath}/spec2zrt'
    proc = await asyncio.create_subprocess_shell(
        f"{code} -d -r {reduceVel} {offset}",
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=runName)

    stdout, stderr = await proc.communicate()

    print(f'[{code} exited with {proc.returncode}]')
    if stdout:
        print(f'[stdout]\n{stdout.decode()}')
    if stderr:
        print(f'[stderr]\n{stderr.decode()}')
    if proc.returncode != 0:
        print(f'Warning, {code} did not exit successfully')
    return proc.returncode

retVal = asyncio.run(mgenkennett(model))
if retVal == 0:
    #asyncio.run(runSpec2zrt(reduceVel, offset))
    results = specfile.readSpecFile(os.path.join(runName, "mspec"),reduceVel = reduceVel, offset = offset)
    for tsObj in results['timeseries']:
        distStr = f"D{tsObj['distance']}"[:5].replace('.','_').strip('_')
        tsObj['depth'] = round(tsObj['depth'], 5);
        commonHeader = {
            'sampling_rate': results['inputs']['frequency']['nyquist']*2.0,
            'channel': 'BHZ',
            'station': distStr,
            'starttime': UTCDateTime(0)+tsObj['timeReduce'],
            'sac': {
                    'b': tsObj['timeReduce'],
                    'dist': tsObj['distance'],
                    'evdp': tsObj['depth']
                }
            }
        # add arrival times and phase name as flags in SAC header
        for idx, a in enumerate(arrivals):
            commonHeader['sac'][f"t{idx}"] = a.time
            commonHeader['sac'][f"kt{idx}"] = a.name
        print(f"Common header: {json.dumps(commonHeader['sac'], indent=4)}")
        header = Stats(commonHeader)
        header.component = 'Z'
        header.npts = len(tsObj['z'])
        z = obspy.Trace(tsObj['z'], header)
        z.write(os.path.join(runName, f"{header.component}_{distStr}_{tsObj['depth']}.sac"), format="SAC")
        header = Stats(commonHeader)
        header.component = 'R'
        header.npts = len(tsObj['r'])
        r = obspy.Trace(tsObj['r'], header)
        r.write(os.path.join(runName, f"{header.component}_{distStr}_{tsObj['depth']}.sac"), format="SAC")
        header = Stats(commonHeader)
        header.component = 'T'
        header.npts = len(tsObj['t'])
        t = obspy.Trace(tsObj['t'], header)
        t.write(os.path.join(runName, f"{header.component}_{distStr}_{tsObj['depth']}.sac"), format="SAC")

print("I'm done.")
