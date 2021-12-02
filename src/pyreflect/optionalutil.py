import sys
import os
import pprint
import json
import math
import tempfile
from io import StringIO
from .earthmodel import EarthModel, list_distances
from .specfile import readSpecFile, AMP_STYLE_VEL, AMP_STYLE_DISP
from .velocitymodel import AK135F, depth_points_from_layers, load_nd_as_depth_points, extend_whole_earth, save_nd
from .stationmetadata import create_fake_metadata, create_stacode_for_dist
from .distaz import DistAz

try:
    import obspy
    import obspy.taup
    import obspy.taup.taup_create
    import obspy.io.sac
    from obspy.core.trace import Stats
    from obspy.core.utcdatetime import UTCDateTime
    obspy_ok = True
except ImportError:
    # obspy not installed
    obspy_ok = False

def check_obspy_import_ok():
    if not obspy_ok:
        raise Error("function requires obspy, but appears not to be installed, http://obspy.org")

ROUND_SLOWNESS_DIGITS = 7

DEPTH_INDEX=3 # index of depth in pierce points output
WAY_BIG=sys.float_info.max
def estimate_for_phases(dist_params, source_depths, phase_list, base_model="ak135", max_depth_offset=200.0):
    """calc travel times to estimate model depth and slowness values"""
    check_obspy_import_ok()
    nd_model_name = base_model
    if nd_model_name == AK135F:
        nd_model_name = 'ak135' # name in obspy
    taumodel = obspy.taup.TauPyModel(model=nd_model_name )
    radiusOfEarth = 6371 # for flat to spherical ray param conversion, should get from model
    maxDepth = 0
    minRayParam = WAY_BIG
    max_red_vel = 0
    maxRayParam = -1
    earliestArrival = None

    for depth_km in source_depths:
        for dist_km in list_distances(dist_params):
            dist_deg = DistAz.kilometersToDegrees(dist_km)
            arrivals = taumodel.get_pierce_points(source_depth_in_km=depth_km,
                                                distance_in_degree=dist_deg,
                                                phase_list=phase_list)
            for a in arrivals:
                if earliestArrival is None or earliestArrival.time > a.time:
                    earliestArrival = a
                minRayParam = min(minRayParam, a.ray_param)
                maxRayParam = max(maxRayParam, a.ray_param)
                for p in a.pierce:
                    maxDepth = max(maxDepth, p[DEPTH_INDEX])
                a_redvel = DistAz.degreesToKilometers(a.distance) / a.time
                if a_redvel > max_red_vel:
                    max_red_vel = a_redvel

    maxDepth = round(math.ceil(maxDepth + max_depth_offset)) # little bit deeper
    minRayParam = minRayParam/radiusOfEarth # need to be flat earth ray params
    maxRayParam = maxRayParam/radiusOfEarth
    if base_model == "ak135" or base_model == AK135F:
        model = EarthModel.loadAk135f(maxDepth)
    elif base_model == "prem":
        model = EarthModel.loadPrem(maxDepth)
    else:
        raise Exception(f"unknown base mode: {base_model}")
    model.distance = dist_params
    model.sourceDepths = source_depths
    model.slowness = {
        "lowcut": round(minRayParam/2, ROUND_SLOWNESS_DIGITS),
        "lowpass": round(minRayParam*0.9, ROUND_SLOWNESS_DIGITS),
        "highpass": round(maxRayParam*1.1, ROUND_SLOWNESS_DIGITS),
        "highcut": round(maxRayParam*1.5, ROUND_SLOWNESS_DIGITS),
        "controlfac": 1.0
        }
    model.extra["earliest_arrival"] = arrival_to_dict(earliestArrival)
    if max_red_vel > 0:
        model.extra["reduce_velocity"] = max_red_vel
    model.extra["phase_list"] = phase_list
    return model

def arrival_to_dict(a):
    return {
        "time": a.time,
        "distance": a.distance,
        "source_depth": a.source_depth,
        "name": a.name,
        "ray_param": a.ray_param,
        "receiver_depth": a.receiver_depth,
        "takeoff_angle": a.takeoff_angle,
        "incident_angle": a.incident_angle,
        "purist_distance": a.purist_distance,
        "purist_name": a.purist_name,
    }

def load_mspec():
    check_obspy_import_ok()
    pass

def arrivals_for_stream(stream):
    check_obspy_import_ok()
    pass

def create_taupymodel(model, extendmodel=AK135F):
    model_name = model.name.split()[0]
    with tempfile.TemporaryDirectory() as output_folder:
        nd_filename = os.path.join(output_folder, model_name + ".nd")
        output_filename = os.path.join(output_folder, model_name + ".npz")
        points = depth_points_from_layers(model.layers)
        extend_points = load_nd_as_depth_points(extendmodel)
        points = extend_whole_earth(points, extend_points)
        save_nd(points, nd_filename)
        save_nd(points, "/Users/crotwell/saved_model.nd")
        mod_create = obspy.taup.taup_create.TauPCreate(input_filename=nd_filename,
                                    output_filename=output_filename)
        mod_create.load_velocity_model()
        tau_model = mod_create.create_tau_model(mod_create.v_mod)
        tau_model.depth_correct(0.0)
        taup = obspy.taup.TauPyModel()
        taup.model = tau_model
    return taup

def mspec_to_stream(rundirectory, model, reduceVel=None, offset=None, phase_list=None, ampStyle=AMP_STYLE_VEL, mspec_filename='mspec'):
    check_obspy_import_ok()
    stream = None
    inv = None
    bandcode = 'B'
    gaincode = 'H'
    loccode = "SY"
    netcode = "XX"
    taupymodel = create_taupymodel(model, extendmodel=AK135F)

    if reduceVel is None and model.extra['reduce_velocity'] is not None:
        reduceVel = model.extra['reduce_velocity']
    elif reduceVel is None:
        reduceVel = 0.0
    if offset is None and model.extra['offset'] is not None:
        offset = model.extra['offset']
    elif offset is None:
        offset = 0.0
    if phase_list is None and model.extra['phase_list'] is not None:
        phase_list = model.extra['phase_list']

    km_to_deg = 180/taupymodel.model.radius_of_planet
    results = readSpecFile(os.path.join(rundirectory, mspec_filename),reduceVel = reduceVel, offset = offset)
    ampStyle = results['inputs']['time']['ampStyle']
    if ampStyle == AMP_STYLE_VEL:
        idep = obspy.io.sac.header.ENUM_VALS['ivel']
    elif ampStyle == AMP_STYLE_DISP:
        idep = obspy.io.sac.header.ENUM_VALS['idisp']
    else:
        idep = obspy.io.sac.header.ENUM_VALS['iunkn']
    for tsObj in results['timeseries']:
        tsObj['depth'] = round(tsObj['depth'], 5)
        stacode = create_stacode_for_dist(tsObj['distance'])
        commonHeader = {
            'sampling_rate': results['inputs']['frequency']['nyquist']*2.0,
            'channel': bandcode+gaincode+'Z',
            'location': loccode,
            'station': stacode,
            'network': netcode,
            'starttime': UTCDateTime(0)+tsObj['timeReduce'],
            'sac': {
                    'b': tsObj['timeReduce'],
                    'dist': tsObj['distance'],
                    'gcarc': tsObj['distance']*km_to_deg,
                    'evdp': tsObj['depth'],
                    'idep': idep
                }
            }
        if tsObj['mech'] != "mij":
            loccode = tsObj['mech'].upper()
        else:
            print(f"mech: {tsObj['mech']}")
            loccode = 'SY'
        commonHeader['location'] = loccode
        if len(phase_list) != 0:
            # add arrival times and phase name as flags in SAC header
            arrivals = taupymodel.get_travel_times(source_depth_in_km=tsObj['depth'],
                                                  distance_in_degree=tsObj['distance']*km_to_deg,
                                                  phase_list=phase_list)
            for idx, a in enumerate(arrivals):
                commonHeader['sac'][f"t{idx}"] = a.time
                commonHeader['sac'][f"kt{idx}"] = a.name
        header = Stats(commonHeader)
        header.component = 'Z'
        header.npts = len(tsObj['z'])
        z = obspy.Trace(tsObj['z'], header)
        header = Stats(commonHeader)
        header.component = 'R'
        header.npts = len(tsObj['r'])
        r = obspy.Trace(tsObj['r'], header)
        header = Stats(commonHeader)
        header.component = 'T'
        header.npts = len(tsObj['t'])
        t = obspy.Trace(tsObj['t'], header)
        if stream is None:
            stream = obspy.Stream(traces=[ z, r, t])
        else:
            stream.append(z)
            stream.append(r)
            stream.append(t)
        metadata = create_fake_metadata(model, loccode, bandcode, gaincode, ampStyle=AMP_STYLE_VEL)
        chan_inv = obspy.read_inventory(StringIO(metadata))
        inv = combine_inventory(inv, chan_inv)

    stream.attach_response(inv)
    return stream, inv

def combine_inventory(inv, to_add):
    if inv is None:
        return to_add
    for add_n in to_add:
        found_n = False
        for inv_n in inv:
            if add_n.code == inv_n.code:
                found_n = True
                found_s = False
                for add_s in add_n:
                    for inv_s in inv_n:
                        if add_s.code == inv_s.code:
                            found_s = True
                            for c in add_s:
                                found_c = False
                                for inv_c in inv_s:
                                    if inv_c.code == c.code and inv_c.location_code == c.location_code:
                                        founc_c = True
                                        break
                                if not found_c:
                                    inv_s.channels.append(c)
                            break
                    if not found_s:
                        inv_n.stations.append(add_s)
                break
        if not found_n:
            inv.append(add_n)
    return inv
