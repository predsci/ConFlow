import os
import sys
from datetime import timedelta
from pathlib import Path
from glob import glob
from scipy.io import FortranFile
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.environ['DEV'], 'sunflower'))
sys.path.append(os.path.join(os.environ['DEV'], 'sunflower/balltracking'))
import balltrack as blt



vpfs = sorted(glob('../tests/vp*.data'))
vtfs = sorted(glob('../tests/vt*.data'))
nfiles = len(vpfs)

# Create list of file series to average
navg = 20
tranges = [[i, i + navg] for i in range(0, nfiles, navg)]

plt.rcParams.update({'font.size': 16})

offset = 0
nx = 1024
ny = 512
for i, trange in enumerate(tranges[offset:]):
    i = i + offset
    dt = i*timedelta(minutes=navg*15)
    vpf_list = vpfs[trange[0]:trange[1]]
    vtf_list = vtfs[trange[0]:trange[1]]

    vps = np.zeros([navg, ny, nx])
    vts = np.zeros([navg, ny, nx])
    for k, (vpf, vtf) in enumerate(zip(vpf_list, vtf_list)):
        vps[k, ...] = np.fromfile(vpf).reshape([ny, nx])
        vts[k, ...] = np.fromfile(vtf).reshape([ny, nx])
    vp = vps.mean(axis=0)
    vt = vts.mean(axis=0)

    lanes = blt.make_lanes(vp, vt, 20, 1)
    plt.figure(figsize=(16, 10))
    plt.imshow(lanes, origin='lower', vmin=0, vmax=8, extent=(0, 360, 180, 0), cmap='gray_r')
    plt.xlabel('Longitude (deg)')
    plt.ylabel('Colatitude (deg)')
    plt.title(f'Field #{i+1} - elapsed time: {str(dt)}')
    plt.tight_layout()
    plt.savefig(os.path.join(os.environ['DATA'], f'ConFlow/sanity_check/lanes_avg_deg_{i+1:06d}.jpg'))
    plt.close('all')

