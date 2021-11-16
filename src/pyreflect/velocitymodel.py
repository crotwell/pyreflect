import pkgutil
import math
try:
    import crustone
    crustone_ok = True
except ImportError:
    # crustone not installed
    crustone_ok = False

def check_crustone_import_ok():
    if not crustone_ok:
        raise Exception("function requires crustone, but appears not to be installed, https://github.com/crotwell/crust-one")

DEFAULT_QP = 1500
DEFAULT_QS = 600
__c1__ = None

def __load_model_nd__(modelname="prem"):
    nd_data = pkgutil.get_data(__name__, f"data/{modelname}.nd")
    if nd_data is None:
        return None
    return load_model_from_nd(nd_data.decode('ascii'))

def load_model_from_nd(ndtext):
    layers = []
    prevDepth = 0
    firstLine = True
    layer = None
    prevLayer = None
    layerType = "crust"
    for line in ndtext.splitlines():
        if line == "mantle" or line == "outer-core" or line == "inner-core":
            layerType = line
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
                "qp": float(qp),
                "qs": float(qs),
                "tp1": 1.0e4,
                "tp2": 0.0001,
                "ts1": 1.0e4,
                "ts2": 0.0001,
                "type": layerType,
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

def createLayer(thick, vp, vs, rho, qp=DEFAULT_QP, qs=DEFAULT_QS, tp1=1.0e4, tp2=0.0001, ts1=1.0e4, ts2=0.0001, type="unknown"):
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
      "ts2": ts2,
      "type": type
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
        "ts2": layer['ts2'],
        "type": layer['type'],
    }

def layersFromAk135f(maxdepth):
    return layersFromModel('ak135fcont', maxdepth)

def layersFromPrem(maxdepth):
    return layersFromModel('prem', maxdepth)

def layersFromModel(modelname, maxdepth):
    premLayers = __load_model_nd__(modelname)
    if premLayers is None:
        with open(modelname, 'r') as nd_file:
            premLayers = nd_file.read()
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

def load_crustone():
    check_crustone_import_ok()
    global __c1__
    if __c1__ is None:
        __c1__ = crustone.parse()
    return __c1__

def modify_crustone(layers, lat, lon):
    check_crustone_import_ok()
    num_crust_layers = 0
    orig_crust_thick = 0
    while layers[num_crust_layers]['type'] == 'crust':
        orig_crust_thick += layers[num_crust_layers]["thick"]
        num_crust_layers += 1
    decapitate = layers[num_crust_layers:]

    c1 = load_crustone()
    profile = c1.find_profile(lat, lon)
    crustone_thick = profile.crust_thick()
    if crustone_thick > orig_crust_thick:
        # shrink topmost mantle layer to accomodate
        # note this takes into account elevation from crustone, effectively
        # increasing the overall model thickness
        # note profile.layers[0].topDepth is negative for high elevation
        extra_thickness = crustone_thick - orig_crust_thick + profile.layers[0].topDepth
        if decapitate[0]["thick"] > extra_thickness:
            decapitate[0]["thick"] = decapitate[0]['thick'] - (extra_thickness)
        else:
            raise Exception(f"top mantle layer is not thick enough to subtract extra crust: orig: {orig_crust_thick} crustone: {crustone_thick}  top mantle: {decapitate[0]['thick']}")
    else:
        decapitate[0]["thick"] = orig_crust_thick - crustone_thick
    merge_layers = []
    for l in profile.layers:
        if l.topDepth == l.botDepth:
            # zero thick layer, skip
            continue
        elif l.botDepth == 6371:
            # halfspace/mantle, skip
            continue
        else:
            nlayer = from_crustone_layer(l)
            merge_layers.append(nlayer)
    for l in decapitate:
        merge_layers.append(l)
    return merge_layers

def from_crustone_layer(layer):
    return createLayer(layer.thick(), layer.vp, layer.vs, layer.rho, type="crust")


def modifyCrustToOceanic(layers):
    """
    modify crust to be more ocean-like
    see https://igppweb.ucsd.edu/~gabi/crust/crust2-averages.txt
    but simplify here to only be 2 layers
    """
    crustLayer1 = layers[0]
    crustLayer2 = layers[1]
    mantleLayer = layers[2]
    if crustLayer1['thick'] != 15.0 and crustLayer2['thick'] != 9.4 and mantleLayer['vp'] != 8.11061:
        raise Error(f"layers do not look like PREM, previously modified? thickness: {crustLayer1['thick']} {crustLayer2['thick']}")
    origCrustThick = crustLayer1['thick'] + crustLayer2['thick']
    layers[0]['thick'] = 6.96
    layers[0]['vp'] = 6.6
    layers[0]['vs'] = 3.65
    layers[0]['rho'] = 2.90
    layers[1]['thick'] = 12.91-6.96
    layers[1]['vp'] = 7.11
    layers[1]['vs'] = 3.91
    layers[1]['rho'] = 3.05
    layers[2]['thick'] = layers[2]['thick'] + origCrustThick - layers[0]['thick'] - layers[1]['thick']
    return layers
