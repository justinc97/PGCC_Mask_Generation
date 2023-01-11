import numpy as np
import pickle
import pandas as pd
import healpy as hp
import pysm3.units as u
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from astropy.table import Table
import multiprocessing
from time import time
import utils 
import os
os.environ[
    "OMP_NUM_THREADS"
] = "64"
np.random.seed(10)

comp = "IQU"
mdir = "/global/project/projectdirs/cmb/data/planck2018/pr3/frequencymaps/"
m, h = hp.read_map(mdir + "HFI_SkyMap_353_2048_R3.01_full.fits",
                                       [c + "_STOKES" for c in comp], dtype = np.float64, h = True)
m <<= u.K_CMB
m = m.to('uK_RJ', equivalencies = u.cmb_equivalencies(353 * u.GHz))
I,Q,U = m
P = Q**2 + U**2
I = I**2

print('map read')
target_url = 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/HFI_PCCS_GCC_R2.02.fits'
PGCC = Table.read(target_url, format = 'fits')
df = PGCC.to_pandas()
df = df[df.FLUX_QUALITY == 1]    # Cut to best quality sources
df = df[df.FLUX_BLENDING == 0]   # Cut out overlapping sources
df = df[df.FLUX_353_CLUMP <= 100]
df = df.reset_index(drop = True) # Reset dataframe indices

output_dir = "/global/cscratch1/sd/justinc/"

def stack_I(l, b, p):
    return(hp.gnomview(I, rot = [l,b,p],
                      reso = 1, xsize = 50,
                      coord = 'G', no_plot = True,
                      return_projected_map = True))
def stack_P(l, b, p):
    return(hp.gnomview(P, rot = [l,b,p],
                      reso = 1, xsize = 50,
                      coord = 'G', no_plot = True,
                      return_projected_map = True))

def repeatI(l,b,p):
    res = []
    for i in range(len(p)):
        res.append(stack_I(l,b,p[i]))
    return np.mean(res, axis = 0)

def repeatP(l,b,p):
    res = []
    for i in range(len(p)):
        res.append(stack_P(l,b,p[i]))
    return np.mean(res, axis = 0)

print('beginning stacking')
t0 = time()
glo = df['GLON']
gla = df['GLAT']
p = np.random.randint(0,360,5)
num_proc = 16
pool = multiprocessing.Pool(processes = num_proc)

process_I = [pool.apply_async(repeatI,
                               args = (i, j, p)) for i, j in zip(glo, gla)]
res_I = [p.get() for p in process_I]

process_P = [pool.apply_async(repeatP,
                               args = (i, j, p)) for i, j in zip(glo, gla)]
res_P = [p.get() for p in process_P]
pool.close()
pool.join()
print(f'Stacking complete in in {(time()-t0)/60}m')
stack = [res_I, res_P]

###------ Uncomment this if you want to save all patches ------###
with open(output_dir + "delete_me_output.pkl", 'wb') as f:
    pickle.dump(stack, f)