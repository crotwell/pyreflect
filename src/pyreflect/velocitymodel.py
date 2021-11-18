import pkgutil
import math
import copy
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
PREM = "prem"
AK135F = 'ak135fcont'

class VelocityModelPoint:
    def __init__(self, depth, vp, vs):
        self.depth = depth
        self.vp = vp
        self.vs = vs
        self.rho = 0
        self.qp = 0
        self.qs = 0
        self.type = "unknown"
    def __str__(self):
        return f"{self.depth} {self.vp} {self.vs} {self.rho}"

class VelocityModelLayer:
    def __init__(self, thick, vp, vs, rho, qp=DEFAULT_QP, qs=DEFAULT_QS, tp1=1.0e4, tp2=0.0001, ts1=1.0e4, ts2=0.0001, type="unknown"):
        self.thick = thick
        self.vp = vp
        self.vp_gradient = 0
        self.vs = vs
        self.vs_gradient = 0
        self.rho = rho
        self.rho_gradient = 0
        self.qp = float(qp)
        self.qs = qs
        self.tp1 = tp1
        self.tp2 = tp2
        self.ts1 = ts1
        self.ts2 = ts2
        self.type = type
    def as_dict(self):
        return {
            "thick": self.thick,
            "vp": self.vp,
            "vp_gradient": self.vp_gradient,
            "vs": self.vs,
            "vs_gradient": self.vs_gradient,
            "rho": self.rho,
            "rho_gradient": self.rho_gradient,
            "qp": self.qp,
            "qs": self.qs,
            "tp1": self.tp1,
            "tp2": self.tp2,
            "ts1": self.ts1,
            "ts2": self.ts2,
            "type": self.type,
        }
    def as_points(self, top_depth):
        top = VelocityModelPoint(top_depth, self.vp, self.vs)
        top.rho = self.rho
        top.qp = self.qp
        top.qs = self.qs
        top.type = self.type
        bot = VelocityModelPoint(top_depth+self.thick, self.vp+self.vp_gradient*self.thick, self.vs+self.vs_gradient*self.thick)
        bot.rho = self.rho+self.rho_gradient*self.thick
        bot.qp = self.qp
        bot.qs = self.qs
        bot.type = self.type
        return top, bot
    @staticmethod
    def from_dict(data):
        v = VelocityModelLayer(data["thick"], data["vp"], data["vs"], data["rho"])
        if "vp_gradient" in data: v.vp_gradient = data["vp_gradient"]
        if "vs_gradient" in data: v.vs_gradient = data["vs_gradient"]
        if "rho_gradient" in data: v.rho_gradient = data["rho_gradient"]
        if "qp" in data: v.qp = data["qp"]
        if "qs" in data: v.qs = data["qs"]
        if "tp1" in data: v.tp1 = data["tp1"]
        if "tp2" in data: v.tp2 = data["tp2"]
        if "ts1" in data: v.ts1 = data["ts1"]
        if "ts2" in data: v.ts2 = data["ts2"]
        if "type" in data: v.type = data["type"]
        return v
    def __copy__(self):
        c = VelocityModelLayer(self.thick, self.vp, self.vs, self.rho, type=self.type)
        c.vp_gradient = self.vp_gradient
        c.vs_gradient = self.vs_gradient
        c.rho_gradient = self.rho_gradient
        c.qp = self.qp
        c.qs = self.qs
        c.tp1 = self.tp1
        c.tp2 = self.tp2
        c.ts1 = self.ts1
        c.ts2 = self.ts2
        return c
    def __str__(self):
        return f"{self.thick} {self.vp} {self.vs} {self.rho}"

def load_nd_as_depth_points(modelname=AK135F):
    nd_data = pkgutil.get_data(__name__, f"data/{modelname}.nd")
    if nd_data is None:
        return None
    ndtext = nd_data.decode('ascii')
    points = []
    layer_type = "crust"
    for line in ndtext.splitlines():
        if line == "mantle" or line == "outer-core" or line == "inner-core":
            layer_type = line
        else:
            depth,vp,vs,rho,qp,qs = line.split()
            p = VelocityModelPoint(float(depth), float(vp), float(vs))
            p.rho = float(rho)
            p.qp = float(qp)
            p.qs = float(qs)
            p.type = layer_type
            points.append(p)
    return points

def layers_from_depth_points(points):
    prev = points[0]
    layers = []
    for point in points[1:]:
        if prev.depth != point.depth:
            thick = point.depth - prev.depth
            vp_gradient = (point.vp-prev.vp)/thick
            vs_gradient = (point.vs-prev.vs)/thick
            rho_gradient = (point.rho-prev.rho)/thick
            layer = VelocityModelLayer(thick, prev.vp, prev.vs, prev.rho, qp=prev.qp, qs=prev.qs, type=prev.type)
            layer.vp_gradient = vp_gradient
            layer.vs_gradient = vs_gradient
            layer.rho_gradient = rho_gradient
            layers.append(layer)
        else:
            # discontinuity
            pass
        prev = point
    return layers
def save_nd(points, filename):
    with open(filename, 'w') as out:
        prev = points[0]
        for p in points:
            if p.type != prev.type:
                out.write(f"{p.type}\n")
            out.write(f"{p.depth} {p.vp} {p.vs} {p.rho}\n")
            prev = p
def extend_whole_earth(points, extend_points):
    out = []
    for p in points:
        out.append(p)
        last = p
    for p in extend_points:
        if p.depth < last.depth:
            # too shallow
            pass
        elif p.depth == last.depth:
            if p.vp != last.vp or p.vs != last.vs:
                # discontinuity, output
                out.append(p)
            else:
                pass
        else:
            out.append(p)
    return out

def depth_points_from_layers(layers):
    points = []
    top_depth = 0
    for l in layers:
        top, bot = l.as_points(top_depth)
        points.append(top)
        if bot.vp != top.vp or bot.vs != top.vs:
            points.append(bot)
        top_depth = bot.depth
    return points

def layersFromEMC(modelname, maxdepth, lat=0, lon=0):
    url = f"http://service.iris.edu/irisws/earth-model/1/line?model={modelname}&lat={lat}&lon={lon}&format=geocsv&nodata=404"
    with urllib.request.urlopen('http://python.org/') as response:
        respOut = response.read()



def layersFromAk135f(maxdepth):
    return layers_from_model(AK135F, maxdepth)

def layersFromPrem(maxdepth):
    return layers_from_model(PREM, maxdepth)

def layers_from_model(modelname, maxdepth):
    points = load_nd_as_depth_points(modelname)
    model_layers = layers_from_depth_points(points)
    depth = 0
    layers = []
    for layer in model_layers:
        botDepth = depth + layer.thick
        if botDepth < maxdepth:
            layers.append(layer)
            depth = botDepth
        else:
            if botDepth > maxdepth:
                splitlayer = copy.deepcopy(layer)
                splitlayer.thick = maxdepth - depth
                layers.append(splitlayer)
            else:
                layers.append(layer)
            bottom_layer = layers[-1]
            halfspace = copy.deepcopy(bottom_layer)
            halfspace.thick = 0.0
            halfspace.vp = bottom_layer.vp+bottom_layer.thick*bottom_layer.vp_gradient
            halfspace.vs = bottom_layer.vs+bottom_layer.thick*bottom_layer.vs_gradient
            halfspace.vp_gradient = 0.0
            halfspace.vs_gradient = 0.0
            layers.append(halfspace)
            break
    return layers

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
    while layers[num_crust_layers].type == 'crust':
        orig_crust_thick += layers[num_crust_layers].thick
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
        if decapitate[0].thick > extra_thickness:
            decapitate[0].thick = decapitate[0].thick - (extra_thickness)
        else:
            raise Exception(f"top mantle layer is not thick enough to subtract extra crust: orig: {orig_crust_thick} crustone: {crustone_thick}  top mantle: {decapitate[0].thick}")
    else:
        decapitate[0].thick = orig_crust_thick - crustone_thick
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
    return VelocityModelLayer(layer.thick(), layer.vp, layer.vs, layer.rho, type="crust")


def modifyCrustToOceanic(layers):
    """
    modify crust to be more ocean-like
    see https://igppweb.ucsd.edu/~gabi/crust/crust2-averages.txt
    but simplify here to only be 2 layers
    """
    crustLayer1 = layers[0]
    crustLayer2 = layers[1]
    mantleLayer = layers[2]
    if crustLayer1.thick != 15.0 and crustLayer2.thick != 9.4 and mantleLayer.vp != 8.11061:
        raise Error(f"layers do not look like PREM, previously modified? thickness: {crustLayer1.thick} {crustLayer2.thick}")
    origCrustThick = crustLayer1.thick + crustLayer2.thick
    layers[0].thick = 6.96
    layers[0].vp = 6.6
    layers[0].vs = 3.65
    layers[0].rho = 2.90
    layers[1].thick = 12.91-6.96
    layers[1].vp = 7.11
    layers[1].vs = 3.91
    layers[1].rho = 3.05
    layers[2].thick = layers[2].thick + origCrustThick - layers[0].thick - layers[1].thick
    return layers
