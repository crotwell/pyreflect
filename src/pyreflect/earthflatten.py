import math
import copy

# Reference radius
R = 6371.0
#
# l is aribtrary constant for EFT, rho_f = rho_s (r/R)^(l+2)
# set l -1 makes the impedence contrast the same for spherical and flattened
l_factor = -1

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

def eft_layer(layer, top_depth, vp_factor=0.1, vs_factor=0.1, R=R, l_factor=l_factor):
    thick = layer.thick
    if thick == 0:
        eft_layer = copy.deepcopy(layer)
        eft_layer.vp = velocity_flat(layer.vp, top_depth)
        eft_layer.vp_gradient = 0
        eft_layer.vs = velocity_flat(layer.vs, top_depth)
        eft_layer.vs_gradient = 0
        eft_layer.rho = density_flat(layer.rho, top_depth, R=R, l_factor=l_factor)
        eft_layer.rho_gradient = 0
        return [ eft_layer ]
    bot_depth = top_depth+thick
    # P wave
    top_vp_s = layer.vp
    bot_vp_s = top_vp_s+layer.vp_gradient*thick
    top_vp_f = velocity_flat(top_vp_s, top_depth, R=R)
    bot_vp_f = velocity_flat(bot_vp_s, top_depth+thick, R=R)
    nnlyrs_vp = math.ceil(abs(bot_vp_f-top_vp_f)/vp_factor)
    # S wave
    top_vs_s = layer.vs
    bot_vs_s = top_vs_s+layer.vs_gradient*thick
    top_vs_f = velocity_flat(top_vs_s, top_depth, R=R)
    bot_vs_f = velocity_flat(bot_vs_s, bot_depth, R=R)
    nnlyrs_vs = math.ceil(abs(bot_vs_f-top_vs_f)/vs_factor)

    nnlyrs = max(nnlyrs_vp, nnlyrs_vs)
    if nnlyrs == 0:
        nnlyrs = 1
    vp = top_vp_f
    delta_vp_f = (bot_vp_f - top_vp_f)/ nnlyrs
    delta_vs_f = (bot_vs_f - top_vs_f)/ nnlyrs
    vs = top_vs_f
    out_layers = []
    prev_depth = depth_flat(top_depth)
    for idx in range(nnlyrs):
        step_vp = delta_vp_f/2 + idx*delta_vp_f
        top_interp_vp_f = top_vp_f + idx*delta_vp_f
        bot_interp_vp_f = top_interp_vp_f + delta_vp_f
        interp_vp_f = (top_interp_vp_f + bot_interp_vp_f)/2.0
        interp_vs_f = top_vs_f + delta_vs_f/2 + idx*delta_vs_f
        eft_layer = copy.deepcopy(layer)
        depth_s = R - radius_for_deltav(top_vp_s, top_depth, bot_vp_s, bot_depth, bot_interp_vp_f-top_vp_f, R=R)
        depth_f = depth_flat(depth_s)
        eft_layer.thick = depth_f - prev_depth
        eft_layer.vp = interp_vp_f
        eft_layer.vp_gradient = 0
        eft_layer.vs = interp_vs_f
        eft_layer.vs_gradient = 0
        eft_layer.rho = density_flat(layer.rho, depth_s, R=R, l_factor=l_factor)
        eft_layer.rho_gradient = 0
        out_layers.append(eft_layer)
        prev_depth = depth_f
    return out_layers

def depth_flat(depth_sph, R=R):
    return R * math.log(R / (R - depth_sph))

def velocity_flat(vel_sph, depth_sph, R=R):
    return R * vel_sph / (R-depth_sph)

def density_flat(density_sph, depth_sph, R=R, l_factor=l_factor):
    return density_sph * math.pow((R-depth_sph)/R, l_factor+2)

def radius_for_deltav(top_v, top_depth, bot_v, bot_depth, delta_v, R=R):
    top_radius = R - top_depth
    bot_radius = R - bot_depth
    # linear factor c, v = top_v + c * (r-top_radius)
    c = (bot_v - top_v ) / (bot_radius - top_radius)
    r = (( top_radius * top_v - c * top_radius * top_radius) /
    (top_radius*delta_v/ R + top_v - c * top_radius ))
    return r
