
from pyreflect import earthmodel, velocitymodel, optionalutil
import obspy
import obspy.taup
import obspy.taup.taup_create

base_model = "ak135"
runName = 'simple'
maxDepth = 800
lat=34.52
lon=92.70


model = earthmodel.EarthModel.loadAk135f(maxDepth)
model.distance = {
    "type": earthmodel.DIST_SINGLE,
    "distance": 1000,
    "azimuth": 45
}
# print the ak135 model as GER style model
print("Unmodified model:")
print(model.name)
print(model.asGER())
print()
print()

# modify model to use crustone profile for (lat,lon). Also takes into account
# crustone elevation to effectively increase radius of earth for extra elevation
# this makes the output model depth range larger. In most cases this probably
# doesn't matter, but might in places like Tibet
tibet = model.crustone(lat, lon)
#tibet = model
# print the modified model
print("Crust1.0 modified model")
print(tibet.name)
print(tibet.asGER())

tibet.export_layers_as_nd("tibet.nd")
tibet.writeToFile("tibet.ger")
