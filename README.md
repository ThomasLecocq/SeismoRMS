# Ground Motion Displacement RMS vs Time

*an example simple tutorial for getting seismic data, computing the power spectral densities, extracting the RMS and plotting*

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3820046.svg)](http://dx.doi.org/10.5281/zenodo.3820046)

Required:

- python
- obspy (and its dependencies)
- pandas
- jupyter
- notebook
- tqdm

this should be easy to set up in a conda env: ``conda create -c conda-forge -n covid python=3.7 obspy pandas jupyter notebook tqdm``

Author: Thomas Lecocq @seismotom, Fred Massin @fmassin

Run it interactively on [mybinder.org](https://mybinder.org/v2/gh/ThomasLecocq/SeismoRMS/master)


## Original Example:
This was the output of the original code shared on Twitter end of March 2020.
The following data shows the effect of the Social Distancing measures from the
Belgian Government (2020-03-16 at midnight, and 2020-03-18 at midday):

![Example image from this code:](img/covid-19_ucc.png)

## Current Example:
The code has evolved and includes several new useful plots for interpreting the
time series.

### Standard ObsPy PPSD plots:

![Standard PPSD:](img/ucc_ppsd.png)

![Standard PPSD_spectrogram:](img/ucc_ppsd_spectrogram.png)

### New plots:

RMS time series:
![RMS_timeseries:](img/BE.UCC..HHZ-4.0-14.0.png)

Changes per day of week & time of day, before and after lockdown:
![RMS_daily_changes:](img/BE.UCC..HHZ-4.0-14.0-daily.png)

Changes per day of week & time of day, before and after lockdown, visualized as
a 24-hour clock:
![RMS_hourly:](img/BE.UCC..HHZ-4.0-14.0-hourly.png)

Colormapped clock plot, each record of the disc is 1 day:
![RMS_clockmap:](img/BE.UCC..HHZ-4.0-14.0-hourmap.png)

Colormapped grid plot, each column of the grid is 1 day:
![RMS_gridmap:](img/BE.UCC..HHZ-4.0-14.0-gridmap.png)
