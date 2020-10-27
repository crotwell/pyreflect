import math

#
# Replaces a layer with many layers approximationg a gradient in a model
#
# Gradient for Rho is .32 * gradient for P, this is probabibly close enough
# for most applications but...
#
# Philip Crotwell, Aug 1995, Perl version
#  based on nawk scripts from George Randall
#
#  tcl version March 1996
#  python version October 2020 (holy cow software lives for a long time)
#
def applyGradient(layers, gradLayerNum, pgrad, sgrad, nlfactor):
#
# layers           list of the layers
# gradLayerNum     layer index to apply the gradient in
# pgrad            Vp gradient to apply
# sgrad            Vs gradient to apply
# nlfactor         approximate layer thickness for the gradient
#
#
    if gradLayerNum > len(layers):
        raise ValueError(f"model doesn't have enough layers, gradLayerNum={gradLayerNum} > model {len(layers)}")
    if gradLayerNum == len(layers)-1:
        raise ValueError(f"can't apply gradient to halfspace, gradLayerNum={gradLayerNum} == model {len(layers)}")
# Revise Model
#
    nnlyr = len(layers)
    ii = 1
    gnumlayers = 0
    layerToReplace = layers[gradLayerNum]
#         how many new layers do we need
    nnl = math.ceil(layerToReplace['thick']/nlfactor)
#
#   here the layer thickness is already smaller than nlfactor so we shouldn't
#   do anything, just set glayer to be layer
#
    if nnl <= 1:
        return layers

    dz = layerToReplace['thick'] / nnl
    dvp = pgrad * layerToReplace['thick'] / nnl
    dvs = sgrad * layerToReplace['thick'] / nnl
# use emprical formula for rho so that rho gradient is .32 of p gradient
    drho = pgrad * .32 * layerToReplace['thick'] / nnl

    roundDigits = 5 # to avoid 3.1229999999 instaed of 3.123

    gradLayers = []
    for i in range(nnl):
        gradLayers.append( {
            "thick": round(dz, roundDigits),
            "vp": round(layerToReplace['vp'] + i * dvp, roundDigits),
            "vs": round(layerToReplace['vs'] + i * dvs, roundDigits),
            "rho": round(layerToReplace['rho'] + i * drho, roundDigits),
            "qp": layerToReplace['qp'],
            "qs": layerToReplace['qs'],
            "tp1": layerToReplace['tp1'],
            "tp2": layerToReplace['tp2'],
            "ts1": layerToReplace['ts1'],
            "ts2": layerToReplace['ts2']
        } )
    preLayers = layers[0:gradLayerNum]
    postLayers = layers[gradLayerNum+1:]
    if len(postLayers) == 1:
        # gradient in layer before halfspace
        #  Set half space parameters equal to prvious layers parameters to avoid
        #  spurious reflections
        botGradLayer = gradLayers[len(gradLayers)-1]
        postLayers = [ {
            "thick": 0.0,
            "vp": botGradLayer['vp'],
            "vs": botGradLayer['vs'],
            "rho": botGradLayer['rho'],
            "qp": botGradLayer['qp'],
            "qs": botGradLayer['qs'],
            "tp1": botGradLayer['tp1'],
            "tp2": botGradLayer['tp2'],
            "ts1": botGradLayer['ts1'],
            "ts2": botGradLayer['ts2']
        } ]

    return preLayers + gradLayers + postLayers
