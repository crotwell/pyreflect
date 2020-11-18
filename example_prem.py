#!/usr/bin/env python3

from pyreflect import earthmodel, distaz, momenttensor

# load model file for PREM down to 85 km depth
model = earthmodel.EarthModel.loadPrem(85)

# pretty print all the model parameters as json,
# this should help if you need to modify things:
model.writeToJsonFile("testprem.json")

# after editing, can load it back in with
othermodel = earthmodel.EarthModel.loadFromJsonFile("testprem.json")
# write out in GER form without the eft, for looking at model
othermodel.writeToFile("modifiedprem.ger")
# and then write back out for actual running synthetics
othermodel.eft().writeToFile("modifiedprem_eft.ger")
# and then run with
# mgenkennett modifiedprem_eft.ger
