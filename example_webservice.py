#!/usr/bin/env python3

from pyreflect import earthmodel, distaz, momenttensor, specfile, stationmetadata
import asyncio
import math
import pprint
import numpy
from io import StringIO, BytesIO

# libcomcat mostly uses same dependencies as obspy, but also needs pyproj

from obspy.clients.fdsn import Client
from obspy.imaging.beachball import MomentTensor
from obspy.taup import TauPyModel
import obspy
import libcomcat
from libcomcat.dataframes import find_nearby_events
from libcomcat.classes import DetailEvent
from libcomcat.search import search
from obspy.io.quakeml.core import Unpickler
from obspy import read
from obspy.core.trace import Stats
from obspy.core.utcdatetime import UTCDateTime

serveSeis = None

try:
    import serveobspy
    serveSeis = serveobspy.ServeObsPy('www')
    serveSeis.serveData()
except ImportError:
    serveSeis = None

# earthquake origin time
#start = obspy.UTCDateTime('2020-10-28T09:02:32')
#start = obspy.UTCDateTime('2020-11-13 09:13:51')
start = obspy.UTCDateTime('2020-11-20T03:19:15')
# eq search paramaters
twindow = 10 # sec
radius = 3 # deg
lat = -53.9
lon = 140.5
minmag = 5

summary_events = search(starttime=start.datetime,
                        endtime = (start+twindow).datetime,
                        maxradius=radius,
                        latitude=lat,
                        longitude=lon,
                        minmagnitude=minmag)
if len(summary_events) == 0:
    print("unable to find event")
    exit(1)

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
else:
    print("unable to fine moment tensor for even")
    exit(1)

if serveSeis is not None:
    serveSeis.catalog = mtcatalog
amp_scale_fac = momenttensor.moment_scale_factor(fm.moment_tensor.scalar_moment)
print(f"mt scale: {fm.moment_tensor.scalar_moment} Nm,  scale fac: {amp_scale_fac}")


searchnetcode = "II"
searchstacode = "TAU"

# find station to calculate distance
sta_client = Client("IRIS")
inventory = sta_client.get_stations(network=searchnetcode, station=searchstacode,
                                level="station",
                                starttime=start,
                                endtime=start+10*60)

network = inventory[0]
station = network.stations[0]

daz = distaz.DistAz(origin.latitude, origin.longitude, station.latitude, station.longitude)

inventory = sta_client.get_stations(network=network.code, station=station.code,
                                location="00", channel="LH?,BH?,HH?",
                                level="response",
                                starttime=start,
                                endtime=start+10*60)
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
print(f"ray param: min {minRayParam}  max {maxRayParam}")

# base model, prem down to 80 km
model = earthmodel.EarthModel.loadPrem(maxDepth)
model.gradientthick = 25 # probably too big
model.eftthick = 50 # probably too big
model.moment_tensor = fm.moment_tensor
model.sourceDepths = [ origin.depth/1000 ] # in km
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

# change to get all ray param
model.slowness['lowcut'] = 0
model.slowness['lowpass'] = 0

#print(model)

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

inventory.write("example.staml", format="stationxml")

outfilebase = "test"
if model.name and len(model.name) > 0:
    outfilebase = model.name.replace(' ', '_')
sphFilename = f"{outfilebase}.ger"
flatFilename = f"{outfilebase}_eft.ger"

# unflattened model (spherical)
model.writeToFile(sphFilename)
# flattened model
model.eft().writeToFile(flatFilename)

# plot beachball with obspy
#from obspy.imaging.beachball import beachball
mtArray = momenttensor.to_beachballarray(fm.moment_tensor.tensor)
print(", ".join(str(x) for x in mtArray))
#fig = beachball(mtArray)

#input("return to quit")

reduceVel = daz.getDistanceKm() / earliestArrival.time
offset = -30



starttime = origin.time + daz.getDistanceKm() / reduceVel + offset
timewidth = model.frequency['numtimepoints'] / (2*model.frequency['nyquist'])

# save waveforms locally
waveforms = sta_client.get_waveforms(network=network.code, station=station.code,
                                location="00", channel="LH?,BH?,HH?",
                                starttime=starttime,
                                endtime=starttime+timewidth)
for tr in waveforms:
    tr.write(f"{tr.id}.sac", format='SAC')
waveforms.attach_response(inventory)

if serveSeis is not None:
    serveSeis.stream = waveforms

async def mgenkennett(model):
    mgenkennett = '../RandallReflectivity/mgenkennett'
    await runCode(mgenkennett, model)


async def runCode(code, model):
    modelFilename = "model_mgen_eft.ger"
    model.eft().writeToFile(modelFilename)
    proc = await asyncio.create_subprocess_shell(
        f"{code} {modelFilename}",
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE)

    stdout, stderr = await proc.communicate()

    print(f'[{code} exited with {proc.returncode}]')
    if stdout:
        print(f'[stdout]\n{stdout.decode()}')
    if stderr:
        print(f'[stderr]\n{stderr.decode()}')


async def runAndPlot(model):
    await mgenkennett(model)
    synthresults = specfile.readSpecFile('mspec', reduceVel=reduceVel, offset=offset)
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
        tr.write(f"{tr.id}.sac", format='SAC')

    both = synthwaveforms+waveforms
    if serveSeis is not None:
        serveSeis.stream = both

    else:
        both.plot()
    print("...switching to interactive mode")
    import code
    code.interact(local=locals())


asyncio.run(runAndPlot(model))
