# 2021-flame-structure-diffusion

This repository contains the data, plotting scripts, and figures associated with the paper "Assessing diffusion model impacts on turbulent transport and flame structure in lean premixed flames" by Aaron J. Fillo, Peter E. Hamlington, and Kyle E. Niemeyer.

These scripts require Matlab (tested on R2021a), but have no external dependencies.

The `plot_*.m` scripts reproduce the figures (`.png` or `.pdf`) using data contained in this archive. The `vorticity_MA.mat` and `vorticity_MC.mat` files are ~2.5 GB and only contained in the Zenodo archive (<https://doi.org/10.5281/zenodo.5146501>).

The `extract_quantities.m` and `calculate_enstrophy_budget.m` scripts read the original NGA
output data files produced in the associated study. This dataset is ~422 GB and archived elsewhere.

The `NGA_grid_reader.m`, `NGAdatareader.m`, and `NGAdatareader_large.m` files were originally 
obtained from members of the [FORCE](https://www.theforce.caltech.edu) research group at Caltech,
and we thank them for their support.

## Attribution

The code in this repository is licensed under the BSD 3-clause license (see the `LICENSE` file
for details), unless otherwise indicated. The figures are licensed under the [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

The Matlab colormaps in the `colormaps` directory were made available by Ander Biguri (2021). Perceptually uniform colormaps (<https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps>), MATLAB Central File Exchange, retrieved July 27, 2021. The original source is <https://bids.github.io/colormap/>.
