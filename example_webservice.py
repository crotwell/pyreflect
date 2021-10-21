#!/usr/bin/env python3

from pyreflect import earthmodel, velocitymodel, distaz, momenttensor, specfile, stationmetadata
import asyncio
import math
import os
import pprint

# we assume obspy, which includes numpy
import numpy
from io import StringIO, BytesIO

import obspy
from obspy.clients.fdsn import Client
from obspy.imaging.beachball import MomentTensor, beachball
from obspy.taup import TauPyModel
from obspy.io.quakeml.core import Unpickler
from obspy import read
from obspy.core.trace import Stats
from obspy.core.utcdatetime import UTCDateTime

# libcomcat mostly uses same dependencies as obspy, but also needs pyproj
# get with: conda install libcomcat
import libcomcat
from libcomcat.dataframes import find_nearby_events
from libcomcat.classes import DetailEvent
from libcomcat.search import search


# name the run, this will create a directory to output lots of files
runName = "tasmania_prem"
if not os.path.exists(runName):
    os.mkdir(runName)

# if using the viewObspy example from seisplotjs:
# see http://crotwell.github.io/seisplotjs/examples/viewobspy/index.html
serveSeis = None

try:
    import serveobspy
    serveSeis = serveobspy.ServeObsPy('www')
    serveSeis.serveData()
except ImportError:
    serveSeis = None

# earthquake origin time
# this event is a 5.5 south of Australia
start = obspy.UTCDateTime('2020-11-20T03:19:15')
# eq search paramaters
twindow = 10 # sec
radius = 3 # deg
lat = -53.9
lon = 140.5
minmag = 5

# use station II.TAU, in Tasmania as it is close
searchnetcode = "II"
searchstacode = "TAU"

# serach for events
summary_events = search(starttime=start.datetime,
                        endtime = (start+twindow).datetime,
                        maxradius=radius,
                        latitude=lat,
                        longitude=lon,
                        minmagnitude=minmag)
if len(summary_events) == 0:
    print("unable to find event")
    exit(1)

# Use the first event (hopefully only) event found
# extract moment tensor and origin parameters
de = summary_events[0].getDetailEvent()
if de.hasProduct('moment-tensor'):
    mtProd = de.getProducts('moment-tensor')
    unpickler = Unpickler()
    qbytes, url = mtProd[0].getContentBytes('xml')
    mtcatalog = unpickler.loads(qbytes)
    event = mtcatalog.events[0]
    origin = event.origins[0]
    fm = event.focal_mechanisms[0]
    #fm = mtcatalog.events[0].preferred_focal_mechanism
    print(f"{origin.time.format_iris_web_service()} {origin.latitude}/{origin.longitude} {fm.moment_tensor.scalar_moment} {fm.moment_tensor.tensor}")
    # save quakeml to file
    mtcatalog.write(os.path.join(runName, "event.qml"), format="QUAKEML")
else:
    print("unable to fine moment tensor for evet")
    exit(1)


# plot beachball with obspy
#from obspy.imaging.beachball import beachball
mtArray = momenttensor.to_beachballarray(fm.moment_tensor.tensor)
print(", ".join(str(x) for x in mtArray))
beachball(mtArray, outfile=os.path.join(runName, "beachball.png"))



if serveSeis is not None:
    serveSeis.catalog = mtcatalog
amp_scale_fac = momenttensor.moment_scale_factor(fm.moment_tensor.scalar_moment)
print(f"mt scale: {fm.moment_tensor.scalar_moment} Nm,  scale fac: {amp_scale_fac}")


# find station to calculate distance
sta_client = Client("IRIS")
inventory = sta_client.get_stations(network=searchnetcode, station=searchstacode,
                                level="station",
                                starttime=start,
                                endtime=start+twindow)
if len(inventory) == 0 or len(inventory[0].stations) == 0:
    raise Error("unable to find station")
network = inventory[0]
station = network.stations[0]

if len(inventory) > 1 or len(inventory[0].stations) > 1:
    print(f"Found more than one net/station, using first: {network.code}.{station.code}")

# calc distance from origin to station
daz = distaz.DistAz(origin.latitude, origin.longitude, station.latitude, station.longitude)

# go get channels with response
inventory = sta_client.get_stations(network=network.code, station=station.code,
                                location="00", channel="LH?,BH?,HH?",
                                level="response",
                                starttime=start,
                                endtime=start+twindow)
inventory.write(os.path.join(runName, f"{network.code}.{station.code}.staxml"), format="STATIONXML")
channels = inventory[0].stations[0].channels

# find LHZ channel to get sample rate
lhz_channel = None
for c in channels:
    if c.code == 'LHZ':
        lhz_channel = c
if lhz_channel is None:
    raise Error(f"Can't find LHZ channel for {network.code}.{station.code}")


if serveSeis is not None:
    serveSeis.inventory=inventory

# calc travel times to estimate model depth and slowness values
taumodel = TauPyModel(model="prem" )
radiusOfEarth = 6371 # for flat to spherical ray param conversion, should get from model

arrivals = taumodel.get_pierce_points(source_depth_in_km=origin.depth/1000,
                                  distance_in_degree=daz.getDistanceDeg(),
                                  phase_list=["P", "S", "SSvmp", "SPvmp", "SSvms"])
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

# might be useful to get wider range of slownesses
minRayParam = 0
# maxRayParam =
print(f"ray param: min {minRayParam}  max {maxRayParam}")


# base model, prem down to 80 km
model = earthmodel.EarthModel.loadPrem(maxDepth)
model.name = runName # optional but useful to set a name
# for the earthquake near Australia set above, and II.TAS, an ocean crust is
# more appropriate. This uses the crust2.0 average oceanic crust just because
model.layers = velocitymodel.modifyCrustToOceanic(model.layers)
model.gradientthick = 25 # probably too big, but quick to run
model.eftthick = 50 # probably too big, but quick to run
model.momentTensor = fm.moment_tensor # use moment tensor from earthquake
model.sourceDepths = [ origin.depth/1000 ] # set source depth in km, quakeml uses m
model.distance = {
    "type": earthmodel.DIST_SINGLE,
    "distance": daz.getDistanceKm(),
    "azimuth": daz.az
}
# this windows around the ray params from the arrivals but note that
# you will get both P and S, so windowing around just an S is not possible
# care with the time window as arrivals wrap due to fft
model.slowness = {
    "lowcut": minRayParam/2,
    "lowpass": minRayParam*0.9,
    "highpass": maxRayParam*1.1,
    "highcut": maxRayParam*2,
    "controlfac": 1.0
    }
model.frequency['numtimepoints'] = 1024
model.frequency['nyquist'] = lhz_channel.sample_rate/2.0
model.frequency['max'] = lhz_channel.sample_rate/2.0

# change to get all ray param, might also change highpass, highcut
model.slowness['lowcut'] = 0
model.slowness['lowpass'] = 0

# model can be printed, but lots of lines of output
#print(model)
# or can save output as json for easier editing
# model.writeToJsonFile("model.json")

bandcode = "L"
gaincode = "H"
loccode = "SY"

synthinventory = stationmetadata.createMetadata(network,
                                                station,
                                                loccode,
                                                bandcode,
                                                gaincode,
                                                model,
                                                fm.moment_tensor.scalar_moment,
                                                specfile.AMP_STYLE_VEL)
print(synthinventory)
synthinventory = obspy.read_inventory(BytesIO(bytes(synthinventory, 'utf-8')))
for c in synthinventory[0].stations[0]:
    print(f"append channel {c}")
    inventory[0].stations[0].channels.append(c)
# force update of inventory as now includes synthethic channels
if serveSeis is not None:
    serveSeis.inventory=inventory

for c in inventory[0].stations[0].channels:
    print(f"channel: {c}")

for c in station.channels:
    print(f"sta channel: {c}")

inventory.write(os.path.join(runName, f"{network.code}.{station.code}_synth.staxml"), format="STATIONXML")

outfilebase = "test"
if model.name and len(model.name) > 0:
    outfilebase = model.name.replace(' ', '_')
sphFilename = f"{outfilebase}.ger"
flatFilename = f"{outfilebase}_eft.ger"

# as json
model.writeToJsonFile(os.path.join(runName, f"{outfilebase}.json"))

# unflattened model (spherical)
model.writeToFile(os.path.join(runName, sphFilename))
# flattened model
model.eft().writeToFile(os.path.join(runName, flatFilename))

#input("return to quit")

reduceVel = daz.getDistanceKm() / earliestArrival.time
offset = -30



starttime = origin.time + daz.getDistanceKm() / reduceVel + offset
timewidth = model.frequency['numtimepoints'] / (2*model.frequency['nyquist'])

# save waveforms locally
print(f" as for waveforms for {network.code}.{station.code}.00.LH?,BH?,HH? {starttime} to {starttime+timewidth}")
waveforms = sta_client.get_waveforms(network=network.code, station=station.code,
                                location="00", channel="LHZ",
                                starttime=starttime,
                                endtime=starttime+timewidth)
for tr in waveforms:
    tr.write(os.path.join(runName, f"{tr.id}.sac"), format='SAC')
waveforms.attach_response(inventory)

if serveSeis is not None:
    serveSeis.stream = waveforms

async def mgenkennett(model):
    mgenkennett = '../../RandallReflectivity/mgenkennett'
    await runCode(mgenkennett, model)


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


async def runAndPlot(model):
    await mgenkennett(model)
    synthresults = specfile.readSpecFile(os.path.join(runName, 'mspec'), reduceVel=reduceVel, offset=offset)
    synthwaveforms = None
    for tsObj in synthresults['timeseries']:
        distStr = f"D{tsObj['distance']}"[:5].replace('.','_').strip('_')
        commonHeader = {
            'sampling_rate': synthresults['inputs']['frequency']['nyquist']*2,
            'network': network.code,
            'station': station.code,
            'channel': bandcode+gaincode+'Z',
            'location': loccode,
            'starttime': UTCDateTime(origin.time)+tsObj['timeReduce'],
            'sac': {
                    'b': tsObj['timeReduce'],
                    'dist': tsObj['distance'],
                    'evdp': tsObj['depth']
                }
            }
        header = Stats(commonHeader)
        header.component = 'Z'
        header.npts = len(tsObj['z'])
        z = obspy.Trace(tsObj['z'], header)
        header = Stats(commonHeader)
        header.component = 'R'
        header.npts = len(tsObj['z'])
        r = obspy.Trace(tsObj['r'], header)
        header = Stats(commonHeader)
        header.component = 'T'
        header.npts = len(tsObj['z'])
        t = obspy.Trace(tsObj['t'], header)
        stream = obspy.Stream(traces=[z, r, t])
        if synthwaveforms is None:
            synthwaveforms = stream
        else:
            synthwaveforms += stream

    print(f"max: {synthwaveforms.max()}   {waveforms.max()}")

    for tr in synthwaveforms:
        tr.write(os.path.join(runName, f"{tr.id}.sac"), format='SAC')

    both = synthwaveforms+waveforms
    if serveSeis is not None:
        serveSeis.stream = both

    else:
        both.plot()
    print("...switching to interactive mode")
    import code
    code.interact(local=locals())


asyncio.run(runAndPlot(model))
