#!/usr/bin/env python3

import os

from pyreflect import earthmodel, velocitymodel, earthflatten


base_model = "ak135"
runName = 'simple'
maxDepth = 500

model = earthmodel.EarthModel.loadAk135f(maxDepth)

(vp, vs, depth) = model.vp_vs_depth()
(grad_vp, grad_vs, grad_depth) = model.evalGradients().vp_vs_depth()
(eft_vp, eft_vs, eft_depth) = model.eft(vp_factor=.05).vp_vs_depth()

import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
plt.ylim(0, 500)
plt.xlim(5, 12)
plt.gca().invert_yaxis()
plt.plot(vp, depth, label='spherical')
plt.plot(grad_vp, grad_depth, label='spherical gradient')
plt.plot(eft_vp, eft_depth, label='eft')
plt.legend()
plt.savefig('vel_model.png')
