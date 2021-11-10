
from pyreflect import earthmodel, velocitymodel


base_model = "ak135"
runName = 'simple'
maxDepth = 100
lat = 34.28
lon = -81.26

model = earthmodel.EarthModel.loadAk135f(maxDepth)
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
c1layers = velocitymodel.modify_crustone(model.layers, lat, lon)
model.layers = c1layers
model.name = model.name+f" modified for Crust 1.0 at {lat}/{lon}, maxDepth: {model.halfspace_depth()}"
# print the modified model
print(model.name)
print(model.asGER())
