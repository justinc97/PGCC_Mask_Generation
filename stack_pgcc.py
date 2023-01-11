#
#
# date: 2022-22-12
# author: Justin Clancy
# Python 3.9
# Copyright (C) 2022 Justin Clancy justinpclancy@gmail.com
#
#
# These are all functions used in the analysis of the accompanying letter for 
# calculating the RMS polarisation fraction of Planck Galactic cold clumps
# using the stacking technique outlined in the letter.
#
#
import numpy as np
mdir_PR4 = "/global/project/projectdirs/cmb/data/planck2020/pla/frequency_maps/Single-frequency/"
mdir_PR3 = "/global/project/projectdirs/cmb/data/planck2018/pr3/frequencymaps/"
PR4      = "HFI_SkyMap_353_2048_R4.00_full.fits"
PR3      = "HFI_SkyMap_353_2048_R3.01_full.fits"

# Read full sky map
import healpy as hp
import pysm3.units as u
comp = "IQU"
m = hp.read_map(mdir_PR4 + PR4, [c + "_STOKES" for c in comp], 
                   dtype = np.float64)
m <<= u.K_CMB
m = m.to('uK_RJ', equivalencies = u.cmb_equivalencies(353 * u.GHz))
m = m ** 2
I,Q,U = m
P = Q + U


import pandas as pd
from astropy.table import Table
target_url = 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/HFI_PCCS_GCC_R2.02.fits'
PGCC = Table.read(target_url, format = 'fits')
df   = PGCC.to_pandas()
df   = df[df.FLUX_QUALITY  == 1]    # Cut to best quality sources
df   = df[df.FLUX_BLENDING == 0]   # Cut out overlapping sources
df   = df.reset_index(drop = True) # Reset dataframe indices

output_dir = "/global/cscratch1/sd/justinc/"


def stack(m,l,b,p):
    return(hp.gnomview(m, rot = [l,b,p],
                  reso = 1, xsize = 50,
                  coord = 'G', no_plot = True,
                  return_projected_map = True))
def repeat(m,l,b):
    p = np.random.randint(0,360,5)
    res = [stack(m,l,b,p[i]) for i in range(len(p))]
    return np.mean(res, axis = 0)

from time import time
import multiprocessing

print("\nStacking unfiltered\n")
t0 = time()
glo = df['GLON']
gla = df['GLAT']
num_proc = 16
pool = multiprocessing.Pool(processes = num_proc)
process_I = [pool.apply_async(repeat, args = (I,i,j)) for i,j in zip(glo,gla)]
process_P = [pool.apply_async(repeat, args = (P,i,j)) for i,j in zip(glo,gla)]
res_I = [p.get() for p in process_I]
res_P = [p.get() for p in process_P]
pool.close()
pool.join()
print(f"Stacking complete in {(time()-t0)/60}m")
stack = [res_I, res_P]
