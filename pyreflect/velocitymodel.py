import pkgutil
import math

def __load_prem__():
    premnd = pkgutil.get_data(__name__, "data/prem.nd").decode('ascii')
    layers = []
    prevDepth = 0
    firstLine = True
    layer = None
    prevLayer = None
    for line in premnd.splitlines():
        if line == "mantle" or line == "outer-core" or line == "inner-core":
            pass
        else:
            depth,vp,vs,rho,qp,qs = line.split()
            depth = float(depth)
            vp = float(vp)
            vs = float(vs)
            rho = float(rho)
            thick = depth-prevDepth
            vp_gradient = 0
            vs_gradient = 0
            if thick > 0:
                vp_gradient = (prevLayer["vp"]-vp)/thick
                vs_gradient = (prevLayer["vs"]-vs)/thick
            layer = {
                "thick": thick,
                "vp": vp,
                "vpgradient": vp_gradient,
                "vs": vs,
                "vsgradient": vs_gradient,
                "rho": rho,
                "qp": qp,
                "qs": qs,
                "tp1": 1.0e4,
                "tp2": 0.0001,
                "ts1": 1.0e4,
                "ts2": 0.0001
            }
            if firstLine:
                firstLine = False
            elif depth == prevDepth:
                pass
            else:
                layers.append(layer)
            prevLayer = layer
            prevDepth = depth
    layers.append(layer)
    return layers

def layersFromEMC(modelname, maxdepth, lat=0, lon=0):
    url = f"http://service.iris.edu/irisws/earth-model/1/line?model={modelname}&lat={lat}&lon={lon}&format=geocsv&nodata=404"
    with urllib.request.urlopen('http://python.org/') as response:
        respOut = response.read()

def createLayer(thick, vp, vs, rho, qp=1000, qs=500, tp1=1.0e4, tp2=0.0001, ts1=1.0e4, ts2=0.0001):
    return {
      "thick": thick,
      "vp": vp,
      "vpgradient": 0.0,
      "vs": vs,
      "vsgradient": 0.0,
      "rho": rho,
      "qp": qp,
      "qs": qs,
      "tp1": tp1,
      "tp2": tp2,
      "ts1": ts1,
      "ts2": ts2
    }

def cloneLayer(layer):
    return {
        "thick": layer["thick"],
        "vp": layer['vp'],
        "vpgradient": layer['vpgradient'],
        "vs": layer['vs'],
        "vsgradient": layer['vsgradient'],
        "rho": layer['rho'],
        "qp": layer['qp'],
        "qs": layer['qs'],
        "tp1": layer['tp1'],
        "tp2": layer['tp2'],
        "ts1": layer['ts1'],
        "ts2": layer['ts2']
    }

def layersFromPrem(maxdepth):
    premLayers = __load_prem__()
    depth = 0
    layers = []
    for layer in premLayers:
        botDepth = depth + layer["thick"]
        if botDepth >= maxdepth:
            if botDepth > maxdepth:
                splitlayer = cloneLayer(layer)
                splitlayer["thick"] = maxdepth - depth
                layers.append(splitlayer)
            else:
                layers.append(layer)
            halfspace = cloneLayer(layer)
            halfspace["thick"] = 0.0
            halfspace["vpgradient"] = 0.0
            halfspace["vsgradient"] = 0.0
            layers.append(halfspace)
            return layers
        else:
            layers.append(layer)
            depth = botDepth
