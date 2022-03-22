from pyreflect import earthmodel, velocitymodel, optionalutil

maxdepth = 100
filename = "copyak135.nd"

points = velocitymodel.load_nd_as_depth_points(filename)
model_layers = velocitymodel.layers_from_depth_points(points)
layers = velocitymodel.trim_layers_for_depth(model_layers, maxdepth)

model = earthmodel.EarthModel()
model.layers = layers
model.name = f"{filename} to {maxdepth}"

print(model.asGER())
