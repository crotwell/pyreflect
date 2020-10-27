import math

roundDigits = 5 # to avoid 3.1229999999 instaed of 3.123

#
# Performs an earth flattening transformation on a model
#
# Philip Crotwell, Aug 1995, Perl version
#  based on nawk scripts from George Randall
#
#  tcl version Feb. 1996
#  python version October 2020
#
def applyEFT(layers, nlfactor):
#
# nlfactor is approximate layer thickness for the EFT
#
# Reference radius
    R = 6371.0
#
# l is aribtrary constant for EFT, rho_f = rho_s (r/R)^(l+2)
# set l -1 makes the impedence contrast the same for spherical and flattened
    l = -1

#
# Revise Model
#
    bots = 0.0
    out = []

    for layerToReplace in layers:
        tops = bots
        bots = bots + layerToReplace['thick']
        topf = R * math.log(R / (R - tops))
        botf = R * math.log(R / (R - bots))
        thickf = botf - topf
        if thickf <= nlfactor:
            nnl = 1
        else:
            nnl = math.ceil(thickf/nlfactor)
        if nnl < 1:
            nnl = 1

        ddz = thickf / nnl
        for ifl in range(nnl):
            layertopf = topf + ifl*ddz
            layerbotf = layertopf + ddz
            expFactor = math.exp((layertopf+layerbotf)/(2*R))
            out.append({
                "thick": round(ddz, roundDigits),
                "vp": round(layerToReplace['vp']*expFactor, roundDigits),
                "vs": round(layerToReplace['vs']*expFactor, roundDigits),
                "rho": round(layerToReplace['rho']*math.pow(expFactor, l+2), roundDigits),
                "qp": layerToReplace['qp'],
                "qs": layerToReplace['qs'],
                "tp1": layerToReplace['tp1'],
                "tp2": layerToReplace['tp2'],
                "ts1": layerToReplace['ts1'],
                "ts2": layerToReplace['ts2']
            })
    # flatten the halfspace
    tops = bots
    topf = R * math.log(R / (R - tops))
    expFactor = math.exp((layertopf+layerbotf)/(2*R))
    out.append({
        "thick": 0.0,
        "vp": round(layerToReplace['vp']*expFactor, roundDigits),
        "vs": round(layerToReplace['vs']*expFactor, roundDigits),
        "rho": round(layerToReplace['rho']*math.pow(expFactor, l+2), roundDigits),
        "qp": layerToReplace['qp'],
        "qs": layerToReplace['qs'],
        "tp1": layerToReplace['tp1'],
        "tp2": layerToReplace['tp2'],
        "ts1": layerToReplace['ts1'],
        "ts2": layerToReplace['ts2']
    })
    return out
