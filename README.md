Polarisation Fraction of \textit{Planck} Galactic Cold Clumps Analysis Code
=============

This is the analysis code and mask generation code created in parallel with the letter to be submitted to MNRAS: Letters in collaboration with the Simons Observatory.

The analysis notebook steps through processes completed in the letter from cutting the \textit{Planck} Galactic Cold Clump Catalogue to the stacking method and error calculations.

The mask generation notebook gives an example of generating a full-sky mask of Galactic cold clumps and applying it to a full-sky frequency map.

Please cite <> if using or referencing any of these results or functions.

## Dependencies
* Python $\ge3.7$
* numpy, pandas, matplotlib
* astropy, healpy
* pysm3 (Python Sky Model 3) -> Can be replaced with just astropy, as this dependency is just for units

## Credits
The models here are built from the data found in the PGCC Catalogue, see arxiv:1502.01599. The code has been built by Justin Clancy with great assistance from Giuseppe Puglisi and completed within the Simons Observatory Collaboration.

## Contents
- `stacking_PGCC_polarisation_analysis_code.ipynb` Works through the letter method and analysis.
- `utils.py` All functions used in the above notebook.
- `generate_masks.ipynb` Works through mask generation and application.
- `generate_mask.py` Functions to generate full-sky Galactic cold clump masks.

If you require use of anything in this work or have any questions, please contact me.
