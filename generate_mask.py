#
#
# date: 2022-16-11
# author: Justin Clancy
# Python 3.9
# Copyright (C) 2022 Justin Clancy justinpclancy@gmail.com
#
#
# These are the functions to generate masks of Galactic cold clumps based on 
# location, size and orientation from the Planck Galactic cold clump catalogue
# PGCC catalogue -> https://irsa.ipac.caltech.edu/data/Planck/html/pgcc_dd.html
#
#

import healpy as hp
import numpy as np
import pandas as pd
from astropy.table import Table
import pysm3.units as u
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import os
os.environ[
    "OMP_NUM_THREADS"
]="64"

def gaussian_source(theta_set, phi_set, a, b, rot):
    """
    Generate a 2D Gaussian profile from a set of HEALPix pixel locations, major
    and minor axes and rotation angle.
    
    -----Inputs-----
    theta_set : angular pixel location
    phi_set   : angular vector location
    a         : FWHM major axis
    b         : FWHM minor axis
    rot       : position/rotation angle
    -----Outputs----
    2D Gaussian elliptical profile
    
    """
    phi_pix, phi_cent = phi_set
    theta_pix, theta_cent = theta_set
    siga = a / np.sqrt(8 * np.log(2))
    sigb = b / np.sqrt(8 * np.log(2))

    A = (phi_pix - phi_cent) * np.cos(rot) - (theta_pix - theta_cent) * np.sin(rot)
    B = (phi_pix - phi_cent) * np.sin(rot) + (theta_pix - theta_cent) * np.cos(rot) 
    return np.exp(-0.5 * ((A**2 / (siga**2)) + (B**2 / (sigb**2))))



def gen_cold_clump_mask(nside = 2048, threshold = 0.9):
    """
    Generate a full-sky mask of Planck Galactic cold clumps from data
    given in the Planck Galactic cold clump catalogue.
    
    -----Inputs-----
    nside     : HEALPix map resolution (default 2048)
    threshold : 2D Gaussian profile zeroing threshold - effectively
                feathers the borders of the source masks, the higher
                the threshold, the sharper the border
    -----Ouputs----
    m         : Full-sky Galactic cold clump mask
    """
    # Load the Planck Galactic cold clump catalogue as a temp file
    target_url = 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/HFI_PCCS_GCC_R2.02.fits'
    # Read catalogue into a dataframe
    PGCC = Table.read(target_url, format = 'fits')
    df = PGCC.to_pandas() 
    """ 
        At this point the function should be edited to make catalogue cuts 
        if needed, i.e., df = df[df['FLUX_QUALITY'] = 1] to restrict to 
        category 1 sources.
    """
    NSources = len(df) # Total number of sources
    npix = hp.nside2npix(nside) # Number of pixels given input resolution
    m = np.ones(npix, dtype = np.float64) # Create a blank full-sky map with value 1
    
    # Read catalogue for size, shape and orientation and define equivalent pixels 
    vecs = hp.ang2vec(df.GLON, df.GLAT, lonlat = True) # Source centre
    fwhm_maj = (np.array(df.GAU_MAJOR_AXIS)*u.arcmin).to_value(u.rad) # Major axis
    fwhm_min = (np.array(df.GAU_MINOR_AXIS)*u.arcmin).to_value(u.rad) # Minor axis
    pos_angle = df.GAU_POSITION_ANGLE # Rotation
    
    # For each source in the catalogue; obtain pixels in a circle of 3 x major axis,
    # turn these pixels into angles of theta and phi along with the sources centre 
    # location, set this region to be equal to the absolute value of 1 minus a 2D
    # Gaussian profile to generate the source mask as a Gaussian centered at zero
    # and rising to 1 at the source extent. Finally, take the source region in the
    # blank map and set to this profile.
    for i in range(NSources):
        pix_region = hp.query_disc(nside = nside,
                                   vec = vecs[i],
                                   radius = 3 * fwhm_maj[i])
        theta_pix, phi_pix = hp.pix2ang(nside, pix_region)
        theta_cent, phi_cent = hp.vec2ang(vecs[i])
        
        # Construct Gaussian elliptical source profile
        profile = abs(1 - gaussian_source((theta_pix, theta_cent),
                                          (phi_pix, phi_cent),
                                          fwhm_maj[i],
                                          fwhm_min[i],
                                          pos_angle[i]))
        for j in range(len(profile)):
            if profile[j] <= threshold:
                profile[j] = 0
        m[pix_region] *= profile
        
    # Return mask
    return m
        