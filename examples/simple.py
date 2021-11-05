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

# path to mgenkennett, relative to simple subdir (may need to add extra ../ )
reflectivityPath = '~/dev/Reflectivity/RandallReflectivity'
distDeg = 45
azimuth = 45
origin_depth_km = 0
reduceVel=8.0  # set below based on earliest arrival
offset=-30
# good idea to have P in phase list to set reduceVel value
phaseList = ["P", "S", "SSvmp", "SPvmp", "SSvms"]
runName = 'simple'
max_freq = 1.0
numPts = 4096


# might be good to scale moment tensor
# assume moment tensor is in newton-meters,
# momenttensor.moment_scale_factor() will convert to
# correct value for GER reflectivity, unity scaled
# to 10**20 dyne-cm
momtensorvals = {
    'm_dd': momenttensor.moment_scale_factor(-0.136e+26),
    'm_nn': momenttensor.moment_scale_factor(-0.073e+26),
    'm_ee': momenttensor.moment_scale_factor(0.209e+26),
    'm_nd': momenttensor.moment_scale_factor(0.074e+26),
    'm_ed': momenttensor.moment_scale_factor(-1.262e+26),
    'm_ne': momenttensor.moment_scale_factor(-0.310e+26)
    }
base_model = "ak135"
#
# End inputs
#
if not os.path.isdir(runName):
    os.mkdir(runName)

# calc travel times to estimate model depth and slowness values
taumodel = TauPyModel(model=base_model )
radiusOfEarth = 6371 # for flat to spherical ray param conversion, should get from model

arrivals = taumodel.get_pierce_points(source_depth_in_km=origin_depth_km,
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


# load model file down to maxDepth km depth

if base_model == "ak135":
    model = earthmodel.EarthModel.loadAk135f(maxDepth)
elif base_model == "prem":
    model = earthmodel.EarthModel.loadPrem(maxDepth)
else:
    raise Exception(f"unknown base mode: {base_model}")
model.name = runName
model.sourceDepths = [origin_depth_km]
model.momentTensor = momtensorvals
model.distance = {
            "type": earthmodel.DIST_SINGLE,
            "distance": distaz.DistAz.degreesToKilometers(distDeg),
            "azimuth": azimuth
        }
model.frequency = {
            "min": 0.0,
            "max": max_freq,
            "nyquist":max_freq,
            "numtimepoints": numPts
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

# runs mgenkennett, note this uses reflectivityPath above to find it
# this is bad, should fix it...
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

# pyreflect can process mspec files natively, so this is only if you want to
# test vs the original fortran spec2zrt
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

# run the synthetic
retVal = asyncio.run(mgenkennett(model))

# if all ok, process the mspec file into sac files
if retVal == 0:
    #asyncio.run(runSpec2zrt(reduceVel, offset))
    results = specfile.readSpecFile(os.path.join(runName, "mspec"),reduceVel = reduceVel, offset = offset)
    for tsObj in results['timeseries']:
        distStr = f"D{tsObj['distance']}"[:6].replace('.','_').strip('_')
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
        arrivals = taumodel.get_pierce_points(source_depth_in_km=tsObj['depth'],
                                              distance_in_degree=tsObj['distance'],
                                              phase_list=phaseList)
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

print(f"{runName}, I'm done.")