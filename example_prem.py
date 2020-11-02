#!/usr/bin/env python3

from pyreflect import earthmodel, distaz, momenttensor

# load model file for PREM down to 85 km depth
model = earthmodel.EarthModel.loadPrem(85)

# write it out, just to look at
model.writeToFile("testprem_nograd.ger")

model = model.evalGradients()

model.writeToFile("testprem.ger")

# apply EFT and write out, suitable for running
model.eft().writeToFile("testprem_eft.ger")
