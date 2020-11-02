#!/usr/bin/env python3

from pyreflect import earthmodel, distaz, momenttensor

# libcomcat mostly uses same dependencies as obspy, but also needs pyproj

from obspy.clients.fdsn import Client
from obspy.imaging.beachball import MomentTensor
import obspy
import libcomcat
from libcomcat.dataframes import find_nearby_events
from libcomcat.classes import DetailEvent
from libcomcat.search import search
from obspy.io.quakeml.core import Unpickler


model = earthmodel.EarthModel.loadPrem(80)

sta_client = Client("IRIS")

start = obspy.UTCDateTime('2020-10-28T09:02:32')
end = start + 20*60

twindow = 10 # sec
radius = 100 # km
lat = -14.9
lon = -75.6
minmag = 5

summary_events = search(starttime=start.datetime,
                        endtime = (start+twindow).datetime,
                        maxradius=radius,
                        latitude=-lat,
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
    print(f"{origin.time.format_iris_web_service()} {origin.latitude}/{origin.longitude}  {fm.moment_tensor.tensor}")
    model.moment_tensor = fm.moment_tensor
    model.sourceDepths.append(origin.depth/1000) # in km
else:
    print("unable to fine moment tensor for even")
    exit(1)


inventory = sta_client.get_stations(network="CO", station="JSC",
                                level="station",
                                starttime=start,
                                endtime=end)

network = inventory[0]
station = network.stations[0]

daz = distaz.DistAz(origin.latitude, origin.longitude, station.latitude, station.longitude)
model.distance = {
    "type": earthmodel.DIST_SINGLE,
    "distance": daz.getDistanceKm(),
    "azimuth": daz.az
}

inventory = sta_client.get_stations(network=network.code, station=station.code,
                                location="00", channel="LH?",
                                level="channel",
                                starttime=start,
                                endtime=end)
channels = inventory[0].stations[0].channels

# unflattened model (spherical)
model.writeToFile("testweb.ger")
# flattened model
model.eft().writeToFile("testweb_eft.ger")

# plot beachball with obspy
#from obspy.imaging.beachball import beachball
mtArray = momenttensor.to_beachballarray(fm.moment_tensor.tensor)
print(", ".join(str(x) for x in mtArray))
#fig = beachball(mtArray)

#input("return to quit")
