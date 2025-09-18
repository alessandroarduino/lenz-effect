# Dynamics of metallic objects in presence of the Lenz effect

This repository collects the data and the scripts used to produce the results collected in the paper "Dynamics of metallic objects in presence of the Lenz effect".

## Overview

The notebook `lenz_effect.ipynb` provides the code necessary to solve the dynamics equations for an aluminum plate moving with one degree of freedom in the static magnetic field of a 1.5 T MRI scanner. The notebook includes the following cases:
   - Circular plate rotating inside the scanner under the action of gravity;
   - Squared plate rotating inside the scanner under the action of gravity;
   - Squared plate translating in the magnet fringe field pushed towards the scanner bore ($y$ = 0);
   - Squared plate translating in the magnet fringe field pushed towards the scanner bore ($y$ = 300 mm).

All the functions defined in the notebook are also listed in the file `lenz_effect.py`, to simplify their adoption in new scripts.
All the functions are thoroughly documented within the files.

An experimental validation of the proposed model and numerical method is provided in the notebook. The script processing the experimental data is contained in folder `experiment`, with the processing outcomes. The experimental data, necessary to run the processing script, can be downloaded from [this Zenodo repository](https://doi.org/10.5281/zenodo.17151932).

## Requirements

The code is provided in the form of Jupyter notebooks, therefore Python and Jupyter are required to operate with them.
Moreover, the following packages are used:
   - [matplotlib](https://matplotlib.org/) (>=3.9.2)
   - [natsort](https://pypi.org/project/natsort/) (>=8.4.0)
   - [numpy](https://numpy.org) (>=2.1.3)
   - [pandas](https://pandas.pydata.org/) (>=2.2.3)
   - [plotly](https://plotly.com/python/) (>=6.3.0)
   - [scikit-image](https://scikit-image.org/) (>=0.25.2)
   - [scipy](https://scipy.org) (>=1.15.2)
   - [seaborn](https://seaborn.pydata.org/) (>=0.13.2)

The package version reported in brackets represents the oldest releases with which the tool has been tested.
Older versions could work as well.

If [xarray](https://docs.xarray.dev/en/stable/index.html) is installed in your environment, be sure that a recent version (e.g., >=2025.9.0) is installed for compatibility with [numpy](https://numpy.org) >=2.0.

## Input data

Two kind of input data are provided and used for the dynamics simulations:
   - The results of electromagnetic simulations, needed to estimate the Lenz effect in rectangular plates for which an analytical solution of the induced currents is not available, are stored in the folder `em_simulations`;
   - The experimental data can be downloaded from [this Zenodo repository](https://doi.org/10.5281/zenodo.17151932). The output of the data processing is already provided in the folder `experiment`.

## Output

The output of the notebook consists of figures that will be generated inside the notebook itself and of tables that will be stored in the folder `results`.
This folder is already populated with the output of the notebook, if you are interested in the data but do not want to generate them again.
