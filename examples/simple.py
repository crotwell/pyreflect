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

from pyreflect import earthmodel, distaz, momenttensor, specfile, optionalutil, velocitymodel

# path to mgenkennett, relative to simple subdir (may need to add extra ../ )
reflectivityPath = '~/dev/Reflectivity/RandallReflectivity'
distDeg = 45
azimuth = 45
origin_depth_km = 0.001 # note cannot have origin on a layer boundary
source_depths = [origin_depth_km]
reduceVel=8.0  # set below based on earliest arrival
offset=-30
# good idea to have P in phase list to set reduceVel value
phase_list = ["P", "S", "SSvmp", "SPvmp", "SSvms"]
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
single_distance = {
            "type": earthmodel.DIST_SINGLE,
            "distance": distaz.DistAz.degreesToKilometers(distDeg),
            "azimuth": azimuth
        }
# or
dist_params = {
    "type": earthmodel.DIST_REGULAR,
    "min": distaz.DistAz.degreesToKilometers(distDeg),
    "delta": distaz.DistAz.degreesToKilometers(2.5),
    "num": 3,
    "azimuth": azimuth
}

# load model file down to maxDepth km depth
source_depths = [origin_depth_km]
model = optionalutil.estimate_for_phases(dist_params, source_depths, phase_list, base_model=velocitymodel.AK135F)
model.name = runName
model.momentTensor = momtensorvals
model.extra['offset'] = -30 # set offset, red vel calculated in model from phases
model.frequency = {
            "min": 0.0,
            "max": max_freq,
            "nyquist":max_freq,
            "numtimepoints": numPts
        }

# pretty print all the model parameters as json,
# this should help if you need to modify things:
model.writeToJsonFile(os.path.join(runName, "simplemodel.json"))
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


# run the synthetic
retVal = asyncio.run(mgenkennett(model))
#retVal = 1

# if all ok, process the mspec file into sac files
if retVal == 0:
    stream, inv = optionalutil.mspec_to_stream(runName, model, ampStyle=specfile.AMP_STYLE_DISP)
    for tr in stream:
        tr.write(os.path.join(runName, tr.id + ".sac"), format="SAC")
    inv.write(os.path.join(runName, "synthetics.stationxml"),
                format="STATIONXML")

print(f"{runName}, I'm done.")
