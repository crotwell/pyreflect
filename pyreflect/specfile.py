import struct
import numpy
import cmath, math

# mspec file is, from fortran read statements:
# fmin, fmax, delf, nfppts, fny, nfpts,nr,nsrc,ndep, azis
# repeat nr ranges
# repeat ndep depths
# repeat ifmin to ifmax: utot, wtot, tntot
#

AMP_STYLE_DISP = 'displacement'
AMP_STYLE_VEL = 'velocity'
SOURCE_STYLE_STEP = 'step'
SOURCE_STYLE_IMPULSE = 'impulse'

MECH_NAMES = {
  0: '_mzz ',
  1: '_mxy ',
  2: '_mxz ',
  3: '_mxx ',
  4: '_myz ',
  5: '_myy '
}


def readSpecFile(filename, reduceVel = 8, offset = -10, ampStyle=AMP_STYLE_VEL, sourceStyle=SOURCE_STYLE_STEP):
    """ Read mspec file generate by mgenkennett or mikjennett

    Output from original Fortran MasterSrcKennett.F:
        write(7) fmin,fmax,delf,nffpts,fny,nfpts,nr,nsrc,nd,azis
        write(7) (rh(i),i=1,nr)
        write(7) (depth(i),i=1,nd)

    """

    results = {}
    inputs = {}
    with open(filename, "rb") as f:
        results['timeseries'] = []
        results['inputs'] = inputs
        inputs['time'] = {
            'reducevel': reduceVel,
            'offset': offset
            }
        struct_fmt = '3fifi3if'
        struct_len = struct.calcsize(struct_fmt)

        startrecord = f.read(4) # fortran starting dummy 4 bytes
        data = f.read(struct_len)
        if not data:
            raise Error(f"tried to read from {filename} but failed")
        endrecord = f.read(4) # fortran endinging dummy 4 bytes
        data_vals = struct.unpack(struct_fmt, data)
        inputs['frequency'] =  {   'min': data_vals[0],
                                    'max': data_vals[1],
                                    'delta': data_vals[2],
                                    'nffpts': data_vals[3],
                                    'nyquist': data_vals[4],
                                    'nfpts': data_vals[5]
                                }
        inputs['numranges'] = data_vals[6]
        inputs['numsources'] = data_vals[7]
        inputs['numdepths'] = data_vals[8]
        inputs['station_azimuth'] = data_vals[9]

        struct_fmt = f"{inputs['numranges']}f"
        struct_len = struct.calcsize(struct_fmt)
        startrecord = f.read(4) # fortran starting dummy 4 bytes
        data = f.read(struct_len)
        endrecord = f.read(4) # fortran endinging dummy 4 bytes
        if not data:
            raise Error(f"tried to read ranges from {filename} but failed")
        inputs['ranges'] = struct.unpack(struct_fmt, data)

        struct_fmt = f"{inputs['numdepths']}f"
        struct_len = struct.calcsize(struct_fmt)
        startrecord = f.read(4) # fortran starting dummy 4 bytes
        data = f.read(struct_len)
        endrecord = f.read(4) # fortran endinging dummy 4 bytes
        if not data:
            raise Error(f"tried to read depths from {filename} but failed")
        inputs['depths'] = struct.unpack(struct_fmt, data)

        freq = inputs['frequency']
        dt = 1. / ( 2 * freq['nyquist'] )
        ifmin = round(freq['min'] / freq['delta'])
        ifmax = round(freq['max'] / freq['delta'])
        ifmax = ifmin+inputs['frequency']['nffpts']-1
        nfppts = ifmax - ifmin + 1
        nft =  2 * ( freq['nfpts'] - 1 )
        for r in range(inputs['numranges']):
            distance = inputs['ranges'][r]
            timeReduce = distance / reduceVel + offset

            for d in range(inputs['numdepths']):
                depth = inputs['depths'][d]
                for s in range(inputs['numsources']):
                    u0 = numpy.zeros(freq['nfpts'], dtype=complex)
                    w0 = numpy.zeros(freq['nfpts'], dtype=complex)
                    tn = numpy.zeros(freq['nfpts'], dtype=complex)

                    for fnum in range(ifmin, ifmax):
                        startrecord = f.read(4) # fortran starting dummy 4 bytes
                        data = f.read(4*2*3)
                        endrecord = f.read(4) # fortran endinging dummy 4 bytes
                        u_real, u_imag, w_real, w_imag, t_real, t_imag = struct.unpack('6f', data)
                        u0[fnum] = complex(u_real,  u_imag)
                        w0[fnum] = complex(w_real,  w_imag)
                        tn[fnum]  = complex(t_real,  t_imag)

                    results['inputs']['u0raw'] = u0.copy()

                    reduceShift = timeReduce * 2*math.pi*freq['delta']
                    

                    for fnum in range(ifmin, ifmax):
                        # apply reducing vel
                        u0[fnum] = u0[fnum] * cmath.exp( complex(0., (fnum)*reduceShift) )
                        w0[fnum] = w0[fnum] * cmath.exp( complex(0., (fnum)*reduceShift) )
                        tn[fnum] = tn[fnum] * cmath.exp( complex(0., (fnum)*reduceShift) )

# Displacement Spectrum has a factor of omeaga**2 from integration of k*dk
#   the Kennett integration is over slowenss p*dp and leaves
#   the remaining factor of omega**2 and and source time spectrum
#   to be included here,
#   Mij moment sources have omega**2 but Fk force sources only omega
#   because the point force source has 1/omega to include
#   Step Source Time Function  is 1/(-i*omega)
#   Impulse Source Time Function is 1
# Velocity spectrum has a factor of (i*omega) * omega**2 [d/dt displ]
                    twopi = 2* math.pi
                    for i in range(freq['nfpts']):
                        fr = (i)*freq['delta']
                        if ampStyle == AMP_STYLE_DISP and sourceStyle == SOURCE_STYLE_STEP:
                      	    ws = 1j * fr * 2 * math.pi
                      	    wsf = 1j
                        elif ampStyle == AMP_STYLE_DISP and sourceStyle == SOURCE_STYLE_IMPULSE:
                      	    ws = -1.*fr*twopi*fr*2 * math.pi
                      	    wsf = -1.*fr*twopi
                        elif ampStyle == AMP_STYLE_VEL and sourceStyle == SOURCE_STYLE_STEP:
                      	    ws = -1.*fr*2 * math.pi*fr*2 * math.pi
                      	    wsf = -1.*fr*2 * math.pi
                        elif ampStyle == AMP_STYLE_VEL and sourceStyle == SOURCE_STYLE_IMPULSE:
                      	    ws = -1j *fr*2 * math.pi*fr*2 * math.pi*fr*2 * math.pi
                      	    wsf = -1j *fr*2 * math.pi*fr*2 * math.pi
                        else:
                            raise Error(f"Dont understand amp/source style: {ampStyle} {sourceStyle}")
                        u0[i] = ws * u0[i]
                        w0[i] = ws * w0[i]
                        tn[i] = ws * tn[i]

                    u0_td = numpy.fft.irfft(u0, nft)
                    w0_td = numpy.fft.irfft(w0)
                    tn_td = numpy.fft.irfft(tn)

                    #scaleFac = 1 /(nft * dt * 4 * math.pi)
                    scaleFac = -1 /( dt * 4 * math.pi) # nft taken care of in fft
                    u0_td = u0_td * scaleFac
                    w0_td = u0_td * scaleFac
                    tn_td = u0_td * scaleFac

                    mech = "mij"
                    greens = {}
                    if (inputs['numsources'] > 1):
                        mech = MECH_NAMES.get(s, "invalid source mech, not 0-5")

                    results['timeseries'].append({
                      "timeReduce": timeReduce,
                      "distance": distance,
                      "depth": depth,
                      "z": u0_td, # z down in GER style synthetics
                      "r": w0_td,
                      "t": tn_td
                    })
    return results
