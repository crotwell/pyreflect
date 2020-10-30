
import json
from .gradient import apply_gradient
from .earthflatten import apply_eft
from .momenttensor import rtp_to_ned

DIST_SINGLE=1
DIST_REGULAR=0
DIST_IRREGULAR=-1

class EarthModel:
    def __init__(self):
        self.isEFT = False
        # layers are:
        # thick vp vs rho qp qs  x x x x
        # I think the last 4 are freq parameters for anisotropy but are not used by the code
        self.layers = [ {
          "thick": 35,
          "vp": 6.5,
          "vs": 3.5,
          "rho": 2.7,
          "qp": 1000,
          "qs": 500,
          "tp1": 1.0e4,
          "tp2": 0.0001,
          "ts1": 1.0e4,
          "ts2": 0.0001
        }, {
          "thick": 0,
          "vp": 8.1,
          "vs": 4.67,
          "rho": 3.32,
          "qp": 2000,
          "qs": 1000,
          "tp1": 1.0e4,
          "tp2": 0.0001,
          "ts1": 1.0e4,
          "ts2": 0.0001
        }]
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
            "type": DIST_SINGLE,
            "distance": 100,
            "azimuth": 45
        }
        # or
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
            "distanceList": [],
            "azimuth": 45
        }
        self.sourceDepths: []
        self.receiverDepth = 0
        self._moment_tensor = {
            "m_nn": 0.0,
            "m_ne": 0.707,
            "m_nd": -0.707,
            "m_ee": 0.0,
            "m_ed": 0.0,
            "m_dd": 0.0
        }
    @staticmethod
    def loadFromFile(filename):
        with open(filename, "r") as f:
            lines = f.readlines()
            return EarthModel.parseGER(lines)
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
                  "vs": float(line[2]),
                  "rho": float(line[3]),
                  "qp": float(line[4]),
                  "qs": float(line[5]),
                  "tp1": 1.0e4,
                  "tp2": 0.0001,
                  "ts1": 1.0e4,
                  "ts2": 0.0001,
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
    def moment_tensor(self):
        return self._moment_tensor
    @setter
    def moment_tensor(self, mt):
        if mt.has('m_rr'):
            self._moment_tensor = rtp_to_ned(mt)
        else:
            self._moment_tensor = mt

    def asJSON(self):
        out = {
            "isEFT": self.isEFT,
            "layers": self.layers,
            "slowness": self.slowness,
            "frequency": self.frequency,
            "distance": self.distance,
            "sourceDepths": self.sourceDepths,
            "receiverDepth": self.receiverDepth,
            "momentTensor": self.momentTensor
        }
        return json.dumps(out, indent=4)
    def asGER(self):
        out = f"{len(self.layers)}\n"
        for l in self.layers:
            out += f"{l['thick']} {l['vp']} {l['vs']} {l['rho']} {l['qp']} {l['qs']} {l['tp1']} {l['tp2']} {l['ts1']} {l['ts2']}\n"
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
        if self.momentTensor:
            out += f"{self.momentTensor['m_nn']} {self.momentTensor['m_ne']} {self.momentTensor['m_nd']} {self.momentTensor['m_ee']} {self.momentTensor['m_ed']} {self.momentTensor['m_dd']}\n"
        return out
    @staticmethod
    def cloneLayer(layer):
        return {
            "thick": layer["thick"],
            "vp": layer['vp'],
            "vs": layer['vs'],
            "rho": layer['rho'],
            "qp": layer['qp'],
            "qs": layer['qs'],
            "tp1": layer['tp1'],
            "tp2": layer['tp2'],
            "ts1": layer['ts1'],
            "ts2": layer['ts2']
        }
    def clone(self):
        out = EarthModel()
        out.isEFT = self.isEFT
        out.layers = [ EarthModel.cloneLayer(x) for x in self.layers ]
        out.slowness = dict(self.slowness)
        out.frequency = dict(self.frequency)
        out.distance = dict(self.distance)
        out.sourceDepths = self.sourceDepths[:]
        out.receiverDepth = self.receiverDepth
        out.momentTensor = dict(self.momentTensor)
        return out
    def gradient(self, gradLayerNum, pgrad, sgrad, nlfactor):
        gradLayers = apply_gradient(self.layers, gradLayerNum, pgrad, sgrad, nlfactor)
        out = self.clone()
        out.layers = gradLayers
        return out
    def eft(self, nlfactor):
        if (self.isEFT):
            raise ValueError("Model has already been flattened")
        eftLayers = apply_eft(self.layers, nlfactor)
        out = self.clone()
        out.isEFT = True
        out.layers = eftLayers
        return out


if __name__ == "__main__":
    with open("../test/testmgen.ger", "r") as f:
        lines = f.readlines()
        mod = EarthModel.parseGER(lines)
        #print(mod.asGER())
        gradmod = mod.gradient(1, .05, .05, 5)
        #print(gradmod.asGER())
        eftmod = gradmod.eft(5)
        print(eftmod.asGER())
        print(eftmod.asJson())
