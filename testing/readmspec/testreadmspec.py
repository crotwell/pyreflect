#!/usr/bin/env python3
import sys
sys.path.append('../../src')

import subprocess
import struct
import pprint
import obspy
import numpy
from obspy.core.trace import Stats
from obspy.core.utcdatetime import UTCDateTime
import matplotlib.pyplot as plt

from pyreflect import earthmodel, distaz, momenttensor, specfile

SPEC2ZRT = "../../../RandallReflectivity/spec2zrt"
# only plot one dist/depth
dist = '1000.0'
depth = '0.001'
reduceVel = 10.0
offset = -60.0

# run spec2zrt
spec2zrt_out = subprocess.run([SPEC2ZRT, '-r', str(reduceVel), str(offset)], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
print(f"Spec2zrt output:")
if len(spec2zrt_out.stderr) > 0:
    raise Exception(f"Probem running spec2zrt: {spec2zrt_out.stderr.decode('utf-8')}")
for line in spec2zrt_out.stdout.decode('utf-8').split('\n'):
    print(f"    {line}")

results = specfile.readSpecFile('mspec',reduceVel = reduceVel, offset = offset)

pp = pprint.PrettyPrinter(width=80, compact=True)
pp.pprint(results['inputs'])

print(f"timeseries: {len(results['timeseries'])}")
traces = []
stream = None

for tsObj in results['timeseries']:
    distStr = f"D{tsObj['distance']}"[:5].replace('.','_').strip('_')
    tsObj['depth'] = round(tsObj['depth'], 5);
    commonHeader = {
        'sampling_rate': results['inputs']['frequency']['nyquist']*2.0,
        'channel': 'BHZ',
        'station': distStr,
        'starttime': UTCDateTime(0)+tsObj['timeReduce'],
        'sac': {
                'b': tsObj['timeReduce'],
                'dist': tsObj['distance'],
                'evdp': tsObj['depth']
            }
        }
    header = Stats(commonHeader)
    header.component = 'Z'
    header.npts = len(tsObj['z'])
    z = obspy.Trace(tsObj['z'], header)
    z.write(f"py_{distStr}_{tsObj['depth']}_{header.component}.sac", format="SAC")
    header = Stats(commonHeader)
    header.component = 'R'
    header.npts = len(tsObj['r'])
    r = obspy.Trace(tsObj['r'], header)
    r.write(f"py_{distStr}_{tsObj['depth']}_{header.component}.sac", format="SAC")
    header = Stats(commonHeader)
    header.component = 'T'
    header.npts = len(tsObj['t'])
    t = obspy.Trace(tsObj['t'], header)
    t.write(f"py_{distStr}_{tsObj['depth']}_{header.component}.sac", format="SAC")
    if tsObj['distance'] == float(dist) and tsObj['depth'] == float(depth):
        stream = obspy.Stream(traces=[z, r, t])
        timeseries_raw = tsObj
    else:
        print(f"Not stream: dist={tsObj['distance']}  depth={tsObj['depth']}")

if stream is None:
    raise Exception(f"Didn't find stream for {dist} {depth}")

#print(stream)
print(stream[0].stats)
#stream.plot()

#exit(0)
tsynth = []
cmp_names = ['z', 'r', 't']
# filename rounds
depth_str = "0.0"
dist_str = "1000"
for cmp_name in cmp_names:
    trace = obspy.read(f"{cmp_name}_{dist_str}_{depth_str}_mij")[0]
    tsynth.append(trace)
    trace.stats.station = 'SAC'
    trace.stats.channel = 'BH'+ cmp_name.upper()
    print(tsynth)
    print(tsynth[0].stats)
#both = stream + tsynth
#both.plot()


for cmp_idx in range(3):
    print(f"##########  {cmp_idx} {cmp_names[cmp_idx].upper()} ###")
    ampS = numpy.abs(numpy.fft.rfft(stream[cmp_idx].data))
    ampT = numpy.abs(numpy.fft.rfft(tsynth[cmp_idx].data))

    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    ax.set_title(f"Amp FFT BH{cmp_names[cmp_idx].upper()}")
    start_point = len(ampS) -10
    start_point = 1
    ax.plot(ampS[start_point:], "r-", label="pyref", linewidth=0.5)
    ax.plot(ampT[start_point:], "b-", label="tsynth", linewidth=0.5)
    plt.legend(loc='lower right')
    plt.savefig(f"myfft_{cmp_idx}.png")
    #plt.show()

    for i in range(len(ampS)-5,len(ampS)):
        print(f"{i} {timeseries_raw['raw']['u0'][i]}   {timeseries_raw['raw']['w0'][i]}  {timeseries_raw['raw']['tn'][i]} ")


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')
    ax.set_title(f"Raw (u0, w0, tn), abs")
    ax.plot(numpy.abs(timeseries_raw['raw']['u0']), "r-", label="u0", linewidth=0.5)
    ax.plot(numpy.abs(timeseries_raw['raw']['w0']), "b-", label="w0", linewidth=0.5)
    ax.plot(numpy.abs(timeseries_raw['raw']['tn']), "g-", label="tn", linewidth=0.5)
    plt.legend(loc='lower right')
    plt.savefig(f"myraw.png")

    # first 100 points:
    only_show_pts_s = 200
    only_show_pts = 400
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(f"BH{cmp_names[cmp_idx].upper()} first {only_show_pts}")
    ax.plot(stream[cmp_idx].times("matplotlib")[:only_show_pts], stream[cmp_idx].data[:only_show_pts], "r-", label="pyref", linewidth=0.5)
    ax.plot(tsynth[cmp_idx].times("matplotlib")[:only_show_pts], tsynth[cmp_idx].data[:only_show_pts], "b-", label="tsynth", linewidth=0.5)
    ax.xaxis_date()
    fig.autofmt_xdate()
    plt.legend(loc='lower right')
    plt.savefig(f"mygraph_{cmp_idx}_first_{only_show_pts}.png")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(f"BH{cmp_names[cmp_idx].upper()}")
    ax.plot(stream[cmp_idx].times("matplotlib"), stream[cmp_idx].data, "r-", label="pyref", linewidth=0.5)
    ax.plot(tsynth[cmp_idx].times("matplotlib"), tsynth[cmp_idx].data, "b-", label="tsynth", linewidth=0.5)
    ax.xaxis_date()
    fig.autofmt_xdate()
    plt.legend(loc='lower right')
    plt.savefig(f"mygraph_{cmp_idx}.png")
    #plt.show()

    max = 0
    maxStream = 0
    maxS_idx = 0
    maxTsynth = 0
    maxT_idx = 0
    maxFFT = 0
    maxFFT_stream = 0
    maxFFT_tsynth = 0
    for i in range(len(ampS)):
        diff = abs(ampS[i]-ampT[i])
        if diff > maxFFT:
            maxFFT = diff
            maxFFT_idx = i
            maxFFT_stream = ampS[i]
            maxFFT_tsynth = ampT[i]
    print(f"Max FFT diff is {maxFFT}, {100.0*maxFFT/maxFFT_tsynth}%, at {maxFFT_idx}, max stream: {maxFFT_stream}, max tsynth: {maxFFT_tsynth}")
    print(f"amp[1024]: {ampS[1024]}  {ampT[1024]}")

    nftp_print = 5
    for i in range(5):
        print(f"{i} {len(ampS)-i-1} {ampS[len(ampS)-i-1]}  {ampT[len(ampT)-i-1]} ")
    print(f"len ampS: {len(ampS)}, len ampT: {len(ampT)}")
    print(f"ampS[{nftp_print}:]: {ampS[(len(ampS)-nftp_print):]}")
    print(f"ampT[{nftp_print}:]: {ampT[(len(ampT)-nftp_print):]}")
    print(f"ampS[{nftp_print}:]: {ampS[len(ampS)-5:]}")
    print(f"ampT[{nftp_print}:]: {ampT[len(ampS)-5:]}")

    if len(stream[0].data) != len(tsynth[cmp_idx].data):
        print(f"Length of data not same: {len(stream[cmp_idx].data)} != {len(tsynth[cmp_idx].data)}")
    for i in range(len(stream[cmp_idx].data)):
        diff = abs(stream[cmp_idx].data[i]-tsynth[cmp_idx].data[i])
        if diff > max:
            max = diff
        if  abs(stream[cmp_idx].data[i]) > abs(maxStream):
            maxStream = stream[cmp_idx].data[i]
            maxS_idx = i
        if  abs(tsynth[cmp_idx].data[i]) > abs(maxTsynth):
            maxTsynth = tsynth[cmp_idx].data[i]
            maxT_idx = i
    print(f"Max diff is {max}, {100.0*max/maxTsynth}%, max stream: {maxStream} ({maxS_idx}), max tsynth: {maxTsynth} ({maxT_idx})")

    diffData = numpy.subtract(stream[cmp_idx].data, tsynth[cmp_idx].data)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(f"Diff BH{cmp_names[cmp_idx].upper()}")
#    ax.plot(stream[0].times("matplotlib")[1:11], diffData[1:11], "r-", label="diff", linewidth=0.5)
    ax.plot(stream[0].times("matplotlib"), diffData, "r-", label="diff", linewidth=0.5)
    ax.xaxis_date()
    fig.autofmt_xdate()
    plt.legend(loc='lower right')
    plt.savefig(f"mydiff_{cmp_idx}.png")
