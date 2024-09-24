# Supporting code for PUF-Bering10K manuscript

This repository holds the supporting code and datasets for the following manuscript:

Kearney, K., Stabeno, P., Hermann, A., Mordy, C. An updated regional model skill assessment for seasonal and interannual variability of bottom temperature across the eastern Bering Sea shelf.  Submitted for review to Frontiers in Marine Science.

## Code

### Matlab analysis code (./popup_analysis_final.m)

The primary script in this repository is the |popup_analysis_final.m| file.  This is a cell-formatted Matlab script that reproduces all post-simulation analysis for this paper, including the generation of all figures in the manuscript.  

This code is included here primarily to support transparency of methods.  However, it should also be runnable by anyone looking to reproduce our results.  A few portions of the workflow involve the use of external software (namely, the [coldpool R package](https://github.com/afsc-gap-products/coldpool) and [ROMSPath offline particle tracking model](https://github.com/imcslatte/ROMSPath)).  We provide the intermediate output of these steps within this repository to allow easier reproduction of our results; flags within the script can be set to either use these pre-generated files or to reproduce them externally.

The code relies on a number of non-native Matlab tools.  Many of these utilities are available in the [MatlabCentral File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/), and can be installed via the Matlab Add-Ons Manager or through direct download and adding to your path:

This author's utilities:

- [labelaxes](https://www.mathworks.com/matlabcentral/fileexchange/171369)
- [mergeaxes](https://www.mathworks.com/matlabcentral/fileexchange/171374)
- [minmax](https://www.mathworks.com/matlabcentral/fileexchange/171379)
- [plotboxpos](https://www.mathworks.com/matlabcentral/fileexchange/9615)
- [plotgrid](https://www.mathworks.com/matlabcentral/fileexchange/171384)
- [shapeprjread](https://www.mathworks.com/matlabcentral/fileexchange/171389)
- [skillstats](https://www.mathworks.com/matlabcentral/fileexchange/171394)

3rd-party utilities:

- [Climate Data Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/70338) (borders, cmocean, ncstruct, shadem, rgb)
- [RunLength](https://www.mathworks.com/matlabcentral/fileexchange/41813)
- [Crameri perceptually uniform scientific colormaps](https://www.mathworks.com/matlabcentral/fileexchange/68546) (crameri)
- [hex2rgb & rgb2hex](https://www.mathworks.com/matlabcentral/fileexchange/46289) (hex2rgb, Note: not required for R2024a+)
- [Hungarian Algorithm for Linear Assignment Problems](https://www.mathworks.com/matlabcentral/fileexchange/20652) (munkres)
- [export_fig](https://www.mathworks.com/matlabcentral/fileexchange/23629)
- [xy2sn, sn2xy](https://www.mathworks.com/matlabcentral/fileexchange/39796)
- [ConsoleProgressBar](https://www.mathworks.com/matlabcentral/fileexchange/30297)
- [arclength](https://www.mathworks.com/matlabcentral/fileexchange/34871)
- [axescoord2figurecoord](https://www.mathworks.com/matlabcentral/fileexchange/13634)
- [distance2curve](https://www.mathworks.com/matlabcentral/fileexchange/34869)
- [interparc](https://www.mathworks.com/matlabcentral/fileexchange/34874)
- [subaxis](https://www.mathworks.com/matlabcentral/fileexchange/3696)

A few additional tools are available through GitHub.  These are part of less polished collections, but the specific functions called in this script are fully documented internally (i.e. via the header text, accessible via the help function).

- [roms](https://github.com/kakearney/roms-pkg) (calcromsz, plotromsrho, stretching)
- [boxmap](https://github.com/kakearney/boxmap-pkg) (boxworldmap)
- [Bering10KPostprocessing](https://github.com/beringnpz/Bering10KPostprocessing) (coldpool_surveyreplicate)

Finally, we use a few functions from the [Gibbs Seawater Toolbox](https://www.teos-10.org/software.htm).

### ROMS standard input (./ROMS_standard_input)

The standard input files for our ROMS Bering10K simulations are also included here.  Again, this is primarily for transparency.  The input files can be paired with our (customized) version of ROMS to reproduce the simulations used in the paper. [The ROMS source code](https://github.com/beringnpz/roms-bering-sea) is available here.  Due to size, we have not included the input forcing files required to run these simulations.  See [Kearney et al., 2020](https://doi.org/10.5194/gmd-13-597-2020) for a detailed description of these inputs. 

## Observational Datasets

The following datasets were used in this paper that are not otherwise available in public archives (yet).

### PMEL Mooring Data (./PMEL_Mooring_Data/)

Mooring data was acquired from PMEL EcoFOCI via Shaun Bell on 12/27/2023.  It includes the mooring data from the Bering Sea, Inner Front, and Kuskokwim deployments, all regularized to a uniform (hourly, ~1m depth-resolution) grid.  Files are netCDF format.  This data is preliminary and is undergoing final quality control in preparation for public archiving.  We provide this version for purposes of reproduction; please contact [TODO](someone's_email?) for an updated version if you wish to use this data for other research.

### Popup Float Data (./PUF_Data)

The popup float data was acquired from PMEL EcoFOCI via Al Hermann on 01/25/2023.  It includes all popup float data from 2021-2022 Bering Sea deployments.  This data is preliminary and is undergoing final quality control in preparation for public archiving.  We provide this version for purposes of reproduction; please contact [TODO](someone's_email?) for an updated version if you wish to use this data for other research.

### NOAA hydrography compendium

Our analysis uses version 0.97 of a multi-decade compendium of NOAA Alaska hydrography:

Pelland, N. A., Nielsen, J. M., Mordy, C. W., Stabeno, P. J., Bell, S. W., Cheng, W., Hermann, A. J., Eisner, L. B., and Gann, J. (2023). The composite Southeast Bering Sea shelf nutricline, within a multi-decade compendium of NOAA Alaska hydrography [Poster presentation]. 2023 Eastern Pacific Ocean Conference, 24-27 September, Fallen Leaf Lake, CA, United States. doi:10.6084/m9.figshare.24556690
(TODO: update citation to Mordy/Pelland manuscript once it is available)

Users wishing to reproduce our calculations can acquire this data from [TODO](final_archive_location) and should adjust the `pellandfol` variable in the `popup_analysis_final.m` script to point to the proper download location.

## Intermediate datasets

### Cold pool raster images (./coldpool_rasters)

Raster images were produced using the [coldpool R package](https://github.com/afsc-gap-products/coldpool) to spatially interpolate survey-sampled model output.  While this calculation can be replicated via that code, we provide the intermediate rasters here for convenience.

### Particle tracking (./particle_tracking)

Input and output files used for the ROMSPath offline particle tracking.