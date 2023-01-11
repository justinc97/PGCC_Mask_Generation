Polarisation Fraction of _Planck_ Galactic Cold Clumps Analysis Code
=============

This is the mask generation code created in parallel with the paper to be submitted to MNRAS in collaboration with the Simons Observatory.

The mask generation notebook gives an example of generating full-sky masks of Galactic cold clumps and applying it to a full-sky frequency map in Stokes I, Q and U.

Please cite <> if using or referencing any of these results or functions.

## Dependencies
* Python $\ge3.7$
* numpy, pandas, matplotlib (just for the notebook)
* astropy, healpy

## Credits
The models here are built from the data found in the PGCC Catalogue, see arxiv:1502.01599. The code has been built by Justin Clancy with great assistance from Giuseppe Puglisi and completed within the Simons Observatory Collaboration.

## Contents
- `generate_masks.ipynb` Works through mask generation and application.
- `generate_mask.py` Functions to generate full-sky Galactic cold clump masks.

If you require use of anything in this work or have any questions, please contact me.
