#
#
# date: 2022-16-11
# author: Justin Clancy
# Python 3.9
# Copyright (C) 2022 Justin Clancy justinpclancy@gmail.com
#
#
# These are all functions used in the analysis of the accompanying letter for 
# calculating the mean squared polarisation fraction of Planck Galactic cold clumps
# using the stacking technique outlined in the letter.
#
#

import numpy as np
import healpy as hp

def sumstack(m):
    """
    Takes a number of sky-patches of equal size and creates a stacked image.
    
    -----Inputs-----
    m          : list of sky-patches
    -----Outputs----
    [SmI, SmP] : Stacked Intensity and Polarisation images
    """
    mI, mP = m
    SmI = np.sum(mI, axis = 0) / len(mI)
    SmP = np.sum(mP, axis = 0) / len(mP)
    return [SmI, SmP]

def BkgAnnulus(data, r, bin_size):
    """
    Creates an annulus of radius, r, and width, bin_size, and 
    returns the mean value within this region. Will also return the
    annulus as a mask.
    
    -----Inputs-----
    data         : Stacked Image
    r            : Radius of annulus (arcmin)
    bin_size     : Width of annulus  (arcmin)
    -----Outputs----
    profile_mean : Mean of pixels within annulus
    mask         : Indices of pixels as a mask
    """
    cx = data.shape[1]//2
    cy = data.shape[0]//2
    a, b = data.shape    
    [X,Y] = np.meshgrid(np.arange(b) - cx, np.arange(a) - cy)
    R = np.sqrt(X**2 + Y**2)
    mask = (np.greater(R, r - bin_size) & np.less(R, r + bin_size))
    values = data[mask]
    test = data.copy()
    test[~mask] *= 0
    profile_mean = np.mean(values)
    profile_std = np.std(values)
    n_pix = len(values)
    return profile_mean, mask

def BkgSubtract(m, r, bin_size):
    """
    Takes intensity and polarisation stacked images and subtracts the
    mean background calculated in an annulus around the central source 
    from the stacked images.
    
    -----Inputs-----
    m          : Intensity and Polarisation stacked sky-patches
    r          : Radius of Annulus (arcmin)
    bin_size   : Width of Annulus (arcmin)
    -----Outputs----
    [I, P] : Stacked Intensity and Polarisation images with background subtracted
    mask   : Returns the annulus indices as a mask for checking if wanted
    """
    I, P = m
    I_annulus_mean, mask = BkgAnnulus(I, r, bin_size) #IS THIS A PROBLEM
    P_annulus_mean, mask = BkgAnnulus(P, r, bin_size)
    I -= I_annulus_mean
    P -= P_annulus_mean
    return [I, P], mask

def bootstrap_stack(sample, sample_size, Bkg_Sub = False):
    """
    Function to take all sky patches and perform bootstrapping with replacement.
    Optional setting to perform background subtraction internally or not.
    
    -----Inputs-----
    sample      : List of all sky patches of equal size (intensity and polarisation)
    sample_size : Number of bootstrap iterations to perform
    Bkg_Sub     : Perform annulus background subtraction after resampled stack 
                  is made (True/False)
    -----Outputs----
    bootstraps  : sample_size number of resulting stacked intensity and polarisation
                  images.
    """
    bootstraps = []
    for _ in range(sample_size):
        src_rawI, src_rawP = sample
        bootstrap_samples = np.random.choice(len(src_rawP), size = len(src_rawP), replace = True)
        bootstrap_P_samples = [src_rawP[i] for i in bootstrap_samples]
        bootstrap_I_samples = [src_rawI[i] for i in bootstrap_samples]
        btstp = [bootstrap_I_samples, bootstrap_P_samples]
        btstp_stacked = sumstack(btstp)     

        if Bkg_Sub == True:
            I, P = BkgSubtract(btstp_stacked, 18, 6)[0]
            bootstrap_PF_samples = [I,P]
           
        else:
            I, P = btstp_stacked
            bootstrap_PF_samples = [I,P]
        bootstraps.append(bootstrap_PF_samples)
    return bootstraps

def sigma_background(data, r, bin_size):
    """
    Calculation of the standard deviation for total squared intensity and polarisation
    computed in the external region of stacked patches.
    
    -----Inputs-----
    data        : Stacked patch
    r           : Radius of background annulus
    bin_size    : Width of background annulus
    -----Outputs----
    profile_std : Standard deviation of external region of stacked patch
    """
    cx = np.array(data).shape[1]//2
    cy = np.array(data).shape[0]//2
    a, b = np.array(data).shape    
    [X,Y] = np.meshgrid(np.arange(b) - cx, np.arange(a) - cy)
    R = np.sqrt(X**2 + Y**2)
    mask = (np.greater(R, r - bin_size) & np.less(R, r + bin_size))
    values = data[mask]
    profile_std = np.std(values)
    return profile_std

def bootstraps2error(bootstraps):
    """
    Takes the bootstrap results of squared intensity and polarisation patches and outputs 
    mean squared results including error bar calculation as per section 3.2 of the associated
    letter.
    
    -----Inputs-----
    bootstraps : sample_size number of bootstrap function outputs
    -----Outputs----
    meanP      : Mean polarisation peak
    meanI      : Mean intesity peak
    sigP       : Standard deviation of polarisation peak
    sigI       : Standard deviation of intensity peak
    B17err     : Error as calculated in letter
    """
    I_stds, P_stds = [], []
    I_cent, P_cent = [], []
    for i in range(len(bootstraps)):
        I_stds.append(sigma_background(bootstraps[i][0], 25, 12))
        P_stds.append(sigma_background(bootstraps[i][1], 25, 12))
        I_cent.append(bootstraps[i][0][25][25])
        P_cent.append(bootstraps[i][1][25][25])
    I_std_mean = np.mean(I_stds)
    P_std_mean = np.mean(P_stds)
    meanP = np.mean(P_cent)
    meanI = np.mean(I_cent)
    sigP  = np.std(P_cent)
    sigI  = np.std(I_cent)
    
    B17err = np.sqrt((meanP/meanI)**2 * ((P_std_mean**2)/(meanP**2) + 
                                         (I_std_mean**2)/(meanI**2)))
    return [meanP, meanI, sigP, sigI, B17err]


def bootstraps2error_bins(bootstrap_bins):
    """
    Takes the bootstrap results of binned squared intensity and polarisation patches and outputs 
    mean squared results including error bar calculation as per section 3.2 of the associated
    letter. This is the same as the 'bootstraps2error' function however it allows for binned data.
    
    -----Inputs-----
    bootstraps : sample_size number of bootstrap function outputs
    -----Outputs----
    meanP      : Mean polarisation peak
    meanI      : Mean intesity peak
    sigP       : Standard deviation of polarisation peak
    sigI       : Standard deviation of intensity peak
    B17err     : Error as calculated in letter
    """
    I_stds, P_stds = [], []
    I_cent, P_cent = [], []
    for i in range(7):
        I_stds_bin, P_stds_bin = [], []
        I_cent_bin, P_cent_bin = [], []
        for j in range(len(bootstrap_bins[0])):
            I_stds_bin.append(sigma_background(bootstrap_bins[i][j][0], 25, 12))
            P_stds_bin.append(sigma_background(bootstrap_bins[i][j][1], 25, 12))
            
            I_cent_bin.append(bootstrap_bins[i][j][0][25][25])
            P_cent_bin.append(bootstrap_bins[i][j][1][25][25])
        I_stds.append(I_stds_bin)
        P_stds.append(P_stds_bin)
        I_cent.append(I_cent_bin)
        P_cent.append(P_cent_bin)
    I_std_mean = np.mean(I_stds, axis = 1)
    P_std_mean = np.mean(P_stds, axis = 1)
    meanP = np.mean(P_cent, axis = 1)
    meanI = np.mean(I_cent, axis = 1)
    sigP  = np.std(P_cent, axis = 1)
    sigI  = np.std(I_cent, axis = 1)
    
    B17err = np.sqrt((meanP/meanI)**2 * ((P_std_mean**2)/(meanP**2) + 
                                         (I_std_mean**2)/(meanI**2)))
    return [meanP, meanI, sigP, sigI, B17err]


def final_results(bootstrapresults):
    """
    Takes the mean squared results from bootstrapping and ouputs the polarisation fraction and error.
    
    -----Inputs-----
    bootstrapresults : Output from 'bootstrap2error' function
    -----Outputs----
    PF               : Mean squared polarisation fraction
    sigPF            : Standard deviation of polarisation fraction
    err              : Error of polarisation fraction as per section 3.2
    """
    P, I, sigP, sigI, err = bootstrapresults
    PF = P/I
    sigPF = sigP/sigI
    return (PF, sigPF, err)

def sqrt_res(squared_results):
    """
    Calculate RMS polarisation fraction and error from mean squared results.
    
    -----Inputs-----
    squared_results : Mean squared polarisation fraction and errors
    -----Outputs----
    np.sqrt(mu)     : RMS polarisation fraction
    sqrt_sig        : RMS standard deviation
    sqrt_err        : RMS polarisation fraction error
    """
    mu, sig, err = squared_results
    sqrt_sig = 0.5 * np.sqrt(mu) * (sig / mu)
    sqrt_err = 0.5 * np.sqrt(mu) * (err / mu)
    return np.sqrt(mu), sqrt_sig, sqrt_err

def meanbinres(binres):
    """
    Calculate the mean squared polarisation fraction and error for binned data.
    
    -----Inputs-----
    binres : Output from binned bootstrapped data
    -----Outputs----
    mean_mu, mean_sig : Total mean polarisation fraction and standard deviation of bins
    mean_muerr, mean_sigerr : Total mean polarisation fraction and error of bins 
    """
    mu, sig, err = binres
    w = sig**(-2)
    w /= np.sum(w)
    werr = err**(-2)
    werr /= np.sum(werr)
    
    mean_mu = np.sum(mu * w)
    mean_muerr = np.sum(mu * werr)
    
    mean_sig = np.sqrt(np.sum(sig**2 * w**2))
    mean_sigerr = np.sqrt(np.sum(err**2 * werr**2))
    return [mean_mu, mean_sig], [mean_muerr, mean_sigerr]