
import pprint
import json
import os
from .gradient import apply_gradient
from .earthflatten import apply_eft
from .momenttensor import rtp_to_ned
from .velocitymodel import layersFromAk135f, layersFromPrem, createLayer, cloneLayer

DIST_SINGLE=1
DIST_REGULAR=0
DIST_IRREGULAR=-1
DEFAULT_MODEL_NAME = "default"

class EarthModel:
    def __init__(self):
        self.name = "default"
        self.gradientthick = 10
        self.eftthick = 5
        self.isEFT = False
        # layers are:
        # thick vp vs rho qp qs  x x x x
        # I think the last 4 are freq parameters for anisotropy but are not used by the code
        self.layers = [ createLayer(35, 6.5, 3.5, 2.7),
                        createLayer(0, 8.1, 4.67, 3.32)]
        self.slowness = {
            "lowcut": 0.005,
            "lowpass": 0.01,
            "highpass": 0.5,
            "highcut": 0.6,
            "controlfac": 1.0
            }
        self.frequency = {
            "min": 0.0,
            "max": 1.0,
            "nyquist":1.0,
            "numtimepoints": 1024
        }
        self.distance = {
            "type": DIST_REGULAR,
            "min": 100,
            "delta": 100,
            "num": 5,
            "azimuth": 45
        }
        #or
        self.distance = {
            "type": DIST_IRREGULAR,
            "distanceList": [100, 150, 330],
            "azimuth": 45
        }
        # or
        self.distance = {
            "type": DIST_SINGLE,
            "distance": 100,
            "azimuth": 45
        }
        self.sourceDepths = [ 0.001 ]
        self.receiverDepth = 0
        self._momentTensor = {
            "m_nn": 0.0,
            "m_ne": 0.707,
            "m_nd": -0.707,
            "m_ee": 0.0,
            "m_ed": 0.0,
            "m_dd": 0.0
        }
    @staticmethod
    def loadPrem(maxdepth ):
        model = EarthModel()
        premLayers = layersFromPrem(maxdepth)
        model.layers = premLayers
        model.name = f"prem to {maxdepth}"
        return model

    @staticmethod
    def loadAk135f(maxdepth ):
        model = EarthModel()
        model.layers = layersFromAk135f(maxdepth)
        model.name = f"ak135f to {maxdepth}"
        return model

    @staticmethod
    def loadFromJsonFile(filename):
        with open(filename, "r") as f:
            jsonModel = json.load(f)
            model = EarthModel.fromDict(jsonModel)
            if model.name == DEFAULT_MODEL_NAME:
                model.name = os.path.basename(filename)
            return model
    def writeToJsonFile(self, filename):
        with open(filename, "w") as f:
            f.write(self.asJSON())
    @staticmethod
    def loadFromFile(filename):
        with open(filename, "r") as f:
            lines = f.readlines()
            model = EarthModel.parseGER(lines)
            model.name = os.path.basename(filename)
            return model
    def writeToFile(self, filename):
        with open(filename, "w") as f:
            f.write(self.asGER())
    @staticmethod
    def parseGER(modelLines):
        out = EarthModel()
        out.layers = []
        i=0
        numLayers = int(modelLines[0].strip())
        for i in range(1, numLayers+1):
            line = modelLines[i].split()
            layer = {
                  "thick": float(line[0]),
                  "vp": float(line[1]),
                  "vpgradient": 0.0,
                  "vs": float(line[2]),
                  "vsgradient": 0.0,
                  "rho": float(line[3]),
                  "qp": float(line[4]),
                  "qs": float(line[5]),
                  "tp1": 1.0e4,
                  "tp2": 0.0001,
                  "ts1": 1.0e4,
                  "ts2": 0.0001
                }
            if len(line) >= 9:
                layer["tp1"]: float(line[6])
                layer["tp2"]: float(line[7])
                layer["ts1"]: float(line[8])
                layer["ts2"]: float(line[9])
            out.layers.append(layer)
        i += 1
        line = modelLines[i].split()
        out.slowness = {
            "lowcut": float(line[0]),
            "lowpass": float(line[1]),
            "highpass": float(line[2]),
            "highcut": float(line[3]),
            "controlfac": float(line[4])
        }
        i += 1
        line = modelLines[i].split()
        out.frequency = {
            "min": float(line[0]),
            "max": float(line[1]),
            "nyquist": float(line[2]),
            "numtimepoints": int(line[3]),
        }
        i += 1
        line = modelLines[i].split()
        disttype = float(line[0])
        azimuth = float(line[1])
        if disttype > 0:
            out.distance = {
                "type": DIST_SINGLE,
                "distance": disttype,
                "azimuth": azimuth
            }
        elif disttype == 0:
            i += 1
            line = modelLines[i].split()
            out.distance = {
                "type": DIST_REGULAR,
                "min": float(line[0]),
                "delta": float(line[1]),
                "num": int(line[2]),
                "azimuth": azimuth
            }
        else: # disttype < 0
            i += 1
            line = modelLines[i].split()
            out.distance = {
                "type": DIST_IRREGULAR,
                "distanceList": [float(x) for x in line],
                "azimuth": azimuth
            }
        i += 1
        line = modelLines[i].split()
        numSources = int(line[0])
        i += 1
        line = modelLines[i].split()
        out.sourceDepths = [ float(x) for x in line ]
        i += 1
        line = modelLines[i].split()
        out.receiverDepth = float(line[0])
        i += 1
        if i < len(modelLines):
            line = modelLines[i].split()
            out.momentTensor = {
                "m_nn": float(line[0]),
                "m_ne": float(line[1]),
                "m_nd": float(line[2]),
                "m_ee": float(line[3]),
                "m_ed": float(line[4]),
                "m_dd": float(line[5])
            }
        return out
    @property
    def momentTensor(self):
        return self._momentTensor
    @momentTensor.setter
    def momentTensor(self, mt):
        tensor = mt
        if 'tensor' in mt:
            tensor = mt.tensor
        if 'm_rr' in tensor:
            self._momentTensor = rtp_to_ned(tensor)
        elif 'm_nd' in tensor:
            self._momentTensor = tensor
        else:
            raise ValueError(f"not sure how to interpret tensor: {tensor}")

    def asDict(self):
        return {
            "name": self.name,
            "gradientthick": self.gradientthick,
            "eftthick": self.eftthick,
            "isEFT": self.isEFT,
            "layers": self.layers,
            "slowness": self.slowness,
            "frequency": self.frequency,
            "distance": self.distance,
            "sourceDepths": self.sourceDepths,
            "receiverDepth": self.receiverDepth,
            "momentTensor": self.momentTensor
        }
    @staticmethod
    def fromDict(data):
        model = EarthModel()
        if "name" in data: model.name = data["name"]
        if "gradientthick" in data: model.gradientthick = data["gradientthick"]
        if "eftthick" in data: model.eftthick = data["eftthick"]
        if "layers" in data: model.layers = data["layers"]
        if "slowness" in data: model.slowness = data["slowness"]
        if "frequency" in data: model.frequency = data["frequency"]
        if "distance" in data: model.distance = data["distance"]
        if "sourceDepths" in data: model.sourceDepths = data["sourceDepths"]
        if "receiverDepth" in data: model.receiverDepth = data["receiverDepth"]
        if "momentTensor" in data: model.momentTensor = data["momentTensor"]
        return model
    def asJSON(self):
        return json.dumps(self.asDict(), indent=4)
    @staticmethod
    def __format_line_as_GER__(items):
        stringItems = [format(x, '.2f') for x in items]
        return " ".join(stringItems)+"\n"

    def asGER(self):
        self.evalGradients()
        out = f"{len(self.layers)}\n"
        for l in self.layers:
            items = [ l['thick'], l['vp'], l['vs'], l['rho'], l['qp'], l['qs'], l['tp1'], l['tp2'], l['ts1'], l['ts2']]
            out += self.__format_line_as_GER__(items)
#            out += f"{format(l['thick'], '.2f')} {l['vp']} {l['vs']} {l['rho']} {l['qp']} {l['qs']} {l['tp1']} {l['tp2']} {l['ts1']} {l['ts2']}\n"
        out += f"{self.slowness['lowcut']} {self.slowness['lowpass']} {self.slowness['highpass']} {self.slowness['highcut']} {self.slowness['controlfac']} \n"
        out += f"{self.frequency['min']} {self.frequency['max']} {self.frequency['nyquist']} {self.frequency['numtimepoints']}\n"
        if self.distance['type'] > 0:
            out += f"{self.distance['distance']} {self.distance['azimuth']}\n"
        elif self.distance['type'] == DIST_REGULAR:
            out += f"{self.distance['type']} {self.distance['azimuth']}\n"
            out += f"{self.distance['min']} {self.distance['delta']} {self.distance['num']}\n"
        else:
            out += f"{self.distance['type']} {self.distance['azimuth']}\n"
            out += f"{len(self.distance['distanceList'])}\n"
            out += f"{' '.join(map(str, self.distance['distanceList']))}\n"
        out += f"{len(self.sourceDepths)}\n"
        out += f"{' '.join(map(str, self.sourceDepths))}\n"
        out += f"{self.receiverDepth}\n"
        if self._momentTensor:
            out += f"{self._momentTensor['m_nn']} {self._momentTensor['m_ne']} {self._momentTensor['m_nd']} {self._momentTensor['m_ee']} {self._momentTensor['m_ed']} {self._momentTensor['m_dd']}\n"
        return out
    def clone(self):
        out = EarthModel()
        out.name = self.name+"Clone"
        out.gradientthick = self.gradientthick
        out.eftthick = self.eftthick
        out.isEFT = self.isEFT
        out.layers = [ cloneLayer(x) for x in self.layers ]
        out.slowness = dict(self.slowness)
        out.frequency = dict(self.frequency)
        out.distance = dict(self.distance)
        out.sourceDepths = self.sourceDepths[:]
        out.receiverDepth = self.receiverDepth
        out.momentTensor = dict(self.momentTensor)
        return out
    def evalGradients(self):
        outLayers = self.layers
        changeMade = True
        while changeMade:
            changeMade = False
            for n in range(len(outLayers)):
                layer = outLayers[n]
                if layer["vpgradient"] != 0.0 or layer["vsgradient"] != 0.0:
                    gradLayers = apply_gradient(outLayers, n, layer['vpgradient'], layer['vsgradient'], self.gradientthick)
                    changeMade = True
                    outLayers = gradLayers
                    break
        out = self.clone()
        out.layers = outLayers
        return out

    def gradient(self, gradLayerNum, pgrad, sgrad, nlfactor):
        gradLayers = apply_gradient(self.layers, gradLayerNum, pgrad, sgrad, nlfactor)
        out = self.clone()
        out.layers = gradLayers
        return out
    def eft(self):
        if (self.isEFT):
            raise ValueError("Model has already been flattened")
        gradMod = self.evalGradients()
        eftLayers = apply_eft(gradMod.layers, self.eftthick)
        out = self.clone()
        out.isEFT = True
        out.layers = eftLayers
        return out
    def __str__(self):
        return pprint.pformat(self.asDict())
