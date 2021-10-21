import math

dyne_cm_per_newton_meter = 1e7 # in dyne cm
randall_unit_scale = 1e20 # in cm per dyne cm
cm_per_m = 100


def moment_scale_factor(scalar_moment_N_m):
    """

c     kennett.f generates regional synthetics following kennett (1983)
c     the units for the kennett synthetic with the specified input
c     units below are:
c     u(cm) = Gij(t)*Fj(t)10-15 + Gij,k*Mjk(t)10-20
c     for F in dynes, M in dyne-cm, * - convolution
c     so a force of 10**15 dynes and a moment of 10**20 dyne-cm
c     are both weighting the respective Greens function or derivative
c     of the Greens function by terms of order 1
c     kennett synthetic will provide displacement in cm with a factor
c     of 10**-20, so Mo of 5*10**20 dyne * cm will be computed
c     by taking 5 times the raw kennett response (for a delta function
c     Mo(t) response) so the 10**-20 scales down the moment figure.
c     (or output in velocity for a step function Mo(t) response)
c     for F(t) in dyne ( or F(w) in dyne * sec )
c     kennett synthetic will provide displacement in cm with a factor
c     of 10**-15

Goal is to return scale factor to take moment and seismogram to m/s amp units
    """
    scale_fac = scalar_moment_N_m * dyne_cm_per_newton_meter / randall_unit_scale / cm_per_m
    return scale_fac

def mw_to_N_m(Mw):
    """
    Mw to Mo conversion from Lay and Wallace p. 384, I assumed that Mo is in
    newton meters hence multiply by 10^7 to change to dyne cm
    (1 Newton = 10^5 dynes and 1 m = 10^2 cm)
    """
    scalar_moment_N_m = math.pow(10,(Mw+10.73)*1.5-7.0)
    return scalar_moment_N_m

def mw_scale_factor(Mw):
    return moment_scale_factor(mw_to_N_m(Mw))

def rtp_to_ned(rtp_momenttensor):
    """Convert a r, theta, phi moment tensor into north, east, down

    r is up
    theta is south
    phi is east

    Due to the way moment tensors work, only the non-diagonal components
    switch sign when reversed.
    """
    mt = rtp_momenttensor

    return {
        "m_dd": mt["m_rr"],
        "m_nn": mt["m_tt"],
        "m_ee": mt["m_pp"],
        "m_nd": mt["m_rt"],
        "m_ed": -1*mt["m_rp"],
        "m_ne": -1*mt["m_tp"]
    }


def ned_to_rtp(ned_momenttensor):
    """Convert a north, east, down moment tensor into r, theta, phi

    r is up
    theta is south
    phi is east

    Due to the way moment tensors work, only the non-diagonal components
    switch sign when reversed.

    USGS, Obspy, Glocal CMT all use rtp
    Kennett and Randall's reflectivity use NED
    """
    mt = rtp_momenttensor
    return {
        "m_rr": mt["m_dd"],
        "m_tt": mt["m_nn"],
        "m_pp": mt["m_ee"],
        "m_rt": mt["m_nd"],
        "m_rp": -1*mt["m_ed"],
        "m_tp": -1*mt["m_ne"]
    }

def to_beachballarray(momenttensor):
    """converts to list suitable for obspy.imaging.beachball
    """
    mt = momenttensor
    if "tensor" in mt:
        mt = mt.tensor
    if "m_nn" in mt:
        mt = ned_to_rtp(mt)
    return [mt.m_rr, mt.m_tt, mt.m_pp, mt.m_rt, mt.m_rp, mt.m_tp]
