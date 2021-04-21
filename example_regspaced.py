#!/usr/bin/env python3

from pyreflect import earthmodel, distaz, momenttensor, stationmetadata

# load model file for PREM down to 85 km depth
model = earthmodel.EarthModel.loadPrem(585)
model.distance = {
    "type": earthmodel.DIST_REGULAR,
    "min": 45*111.19,
    "delta": 5*111.19,
    "num": 5,
    "azimuth": 45
}
loccode="SY"
bandcode = "L"
gaincode = "H"
scalar_moment_N_m = 1e22
ampStyle = stationmetadata.AMP_STYLE_DISP
fakeXml = stationmetadata.createFakeMetadata(model, loccode, bandcode, gaincode, scalar_moment_N_m, ampStyle)

print(fakeXml)

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
