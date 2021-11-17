#!/usr/bin/env python3

import os
import copy

from pyreflect import earthmodel, velocitymodel, earthflatten


base_model = "ak135"
runName = 'simple'
maxDepth = 800
lat=34.52
lon=92.70

model = earthmodel.EarthModel.loadAk135f(maxDepth)
tibet = model.crustone(lat, lon)

points = velocitymodel.load_nd_as_depth_points("ak135fcont")
p_depths = []
p_vp = []
for p in points:
    if p.depth <= maxDepth:
        p_depths.append(p.depth)
        p_vp.append(p.vp)

(vp, vs, depth) = model.vp_vs_depth()
#(grad_vp, grad_vs, grad_depth) = model.evalGradients().vp_vs_depth()
(eft_vp, eft_vs, eft_depth) = model.eft(vp_factor=.05).vp_vs_depth()

tibet_vp, tibet_vs, tibet_depth = tibet.vp_vs_depth()
(eft_tibet_vp, eft_tibet_vs, eft_tibet_depth) = tibet.eft(vp_factor=.05).vp_vs_depth()

eft_p_depths = []
eft_p_vp = []
for p in points:
    if p.depth <= maxDepth:
        eft_p_depths.append(earthflatten.depth_flat(p.depth))
        eft_p_vp.append(earthflatten.velocity_flat(p.vp, p.depth))

import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
#plt.ylim(0, 800)
#plt.xlim(10, 12)
plt.scatter(p_vp, p_depths, label='ak135')
plt.scatter(eft_p_vp, eft_p_depths, label='eft')
plt.gca().invert_yaxis()
plt.plot(vp, depth, label='ak135')
plt.plot(eft_vp, eft_depth, label='eft')
plt.plot(tibet_vp, tibet_depth, label='Tibet')
plt.plot(eft_tibet_vp, eft_tibet_depth, label='eft Tibet')
plt.legend()
plt.savefig('vel_model.png')
