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

DEF_REDUCE_VEL=8.0
DEF_OFFSET=0.0

MECH_NAMES = {
  0: 'zz',
  1: 'xy',
  2: 'xz',
  3: 'xx',
  4: 'yz',
  5: 'yy'
}

def load_specfile(filename):
    """ Read mspec file generate by mgenkennett or mikjennett

    Output from original Fortran MasterSrcKennett.F:
        write(7) fmin,fmax,delf,nffpts,fny,nfpts,nr,nsrc,nd,azis
        write(7) (rh(i),i=1,nr)
        write(7) (depth(i),i=1,nd)

    """

    results = {}
    inputs = {}
    results['inputs'] = inputs
    results['timeseries'] = []

    struct_fmt = '3fifi3if'
    struct_len = struct.calcsize(struct_fmt)

    with open(filename, "rb") as f:
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
        ifmax = ifmin+inputs['frequency']['nffpts']-1
        nfppts = ifmax - ifmin + 1
        nft =  2 * ( freq['nfpts'] - 1 )
        for r in range(inputs['numranges']):
            distance = inputs['ranges'][r]

            for d in range(inputs['numdepths']):
                depth = inputs['depths'][d]
                for s in range(inputs['numsources']):
                    u0 = numpy.zeros(freq['nfpts'], dtype=complex)
                    w0 = numpy.zeros(freq['nfpts'], dtype=complex)
                    tn = numpy.zeros(freq['nfpts'], dtype=complex)

                    for fnum in range(ifmin, ifmax+1):
                        startrecord = f.read(4) # fortran starting dummy 4 bytes
                        data = f.read(4*2*3)
                        endrecord = f.read(4) # fortran endinging dummy 4 bytes
                        u_real, u_imag, w_real, w_imag, t_real, t_imag = struct.unpack('6f', data)
                        u0[fnum] = complex(u_real,  u_imag)
                        w0[fnum] = complex(w_real,  w_imag)
                        tn[fnum]  = complex(t_real,  t_imag)

                    mech = "mij"
                    greens = {}
                    if (inputs['numsources'] > 1):
                        mech = MECH_NAMES.get(s, "invalid source mech, not 0-5")
                    timeseries = {
                      "timeReduce": None,
                      "distance": distance,
                      "depth": depth,
                      "mech": mech,
                      "z": None, # z down in GER style synthetics
                      "r": None,
                      "t": None,
                      "raw": {
                        "u0": u0,
                        "w0": w0,
                        "tn": tn
                      }
                    };
                    results['timeseries'].append(timeseries)
    return results

def to_time_domain(results, reduceVel = DEF_REDUCE_VEL, offset = DEF_OFFSET, ampStyle=AMP_STYLE_VEL, sourceStyle=SOURCE_STYLE_STEP):
    if reduceVel is None:
        reduceVel = DEF_REDUCE_VEL
    if offset is None:
        offset = DEF_OFFSET
    inputs = results['inputs']
    inputs['time'] = {
        'reducevel': reduceVel,
        'offset': offset,
        'ampStyle': ampStyle,
        'sourceStyle': sourceStyle
        }

    freq = inputs['frequency']
    dt = 1. / ( 2 * freq['nyquist'] )
    ifmin = round(freq['min'] / freq['delta'])
    ifmax = ifmin+inputs['frequency']['nffpts']-1
    nft =  2 * ( freq['nfpts'] - 1 )
    for ts in results['timeseries']:
        u0 = ts['raw']['u0'].copy()
        w0 = ts['raw']['w0'].copy()
        tn = ts['raw']['tn'].copy()

        distance = ts['distance']
        timeReduce = distance / reduceVel + offset
        reduceShift = timeReduce * 2*math.pi*freq['delta']

        ts["timeReduce"] = timeReduce

        for fnum in range(ifmin, ifmax+1):
            # apply reducing vel
            u0[fnum] *= cmath.exp( complex(0., (fnum)*reduceShift) )
            w0[fnum] *= cmath.exp( complex(0., (fnum)*reduceShift) )
            tn[fnum] *= cmath.exp( complex(0., (fnum)*reduceShift) )

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
            u0[i] *= ws
            w0[i] *= ws
            tn[i] *= ws

            u0_td = numpy.fft.irfft(u0, nft)
            w0_td = numpy.fft.irfft(w0, nft)
            tn_td = numpy.fft.irfft(tn, nft)

            #scaleFac = 1 /(nft * dt * 4 * math.pi)
            scaleFac = -1 /( dt * 4 * math.pi) # nft taken care of in fft
            ts['z'] = u0_td * scaleFac
            ts['r'] = w0_td * scaleFac
            ts['t'] = tn_td * scaleFac
    return results

def readSpecFile(filename, reduceVel = DEF_REDUCE_VEL, offset = -10.0, ampStyle=AMP_STYLE_VEL, sourceStyle=SOURCE_STYLE_STEP):
    results = load_specfile(filename)
    return to_time_domain(results, reduceVel=reduceVel, offset=offset,ampStyle=ampStyle,sourceStyle=sourceStyle)
