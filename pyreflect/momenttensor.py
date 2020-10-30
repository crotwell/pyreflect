

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
      m_dd: mt.m_rr,
      m_nn: mt.m_tt,
      m_ee: mt.m_pp,
      m_nd: mt.m_rt,
      m_ed: -1*mt.m_rp,
      m_ne: -1*mt.m_tp
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
      m_rr: mt.m_dd,
      m_tt: mt.m_nn,
      m_pp: mt.m_ee,
      m_rt: mt.m_nd,
      m_rp: -1*mt.m_ed,
      m_tp: -1*mt.m_ne
    }
