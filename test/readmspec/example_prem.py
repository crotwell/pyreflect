#!/usr/bin/env python3
import sys
sys.path.append('../..')
from pyreflect import earthmodel, distaz, momenttensor

# load model file for PREM down to 85 km depth
model = earthmodel.EarthModel.loadAk135f(585)
print(f"dist: {model.distance}")
model.distance['distance'] = 1000.0
model.slowness['highpass'] = 0.7
model.slowness['highcut'] = 0.75
mo_N_m = momenttensor.mw_to_N_m(5.93)
rtp_momenttensor = {
"m_rr":-859500000000000000/mo_N_m,
"m_tt":-196500000000000000/mo_N_m,
"m_pp": 1056000000000000000/mo_N_m,
"m_rt": 107300000000000000/mo_N_m,
"m_rp": 55400000000000000/mo_N_m,
"m_tp": 137700000000000000/mo_N_m
}
model.momentTensor = momenttensor.rtp_to_ned(rtp_momenttensor)
#model.momentTensor = {
#"m_nn": 0.0,
#"m_ne": 0.707,
#"m_nd": -0.707,
#"m_ee": 0.0,
#"m_ed": 0.0,
#"m_dd": 0.0
#}
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
