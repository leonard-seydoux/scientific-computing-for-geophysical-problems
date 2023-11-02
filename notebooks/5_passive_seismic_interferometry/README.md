# 5. Cross-correlation of background seismic noise for Green's function retrieval

<img width=150px src="https://upload.wikimedia.org/wikipedia/fr/thumb/1/16/Institut_de_physique_du_globe_de_paris_logo.svg/1200px-Institut_de_physique_du_globe_de_paris_logo.svg.png" />

> This notebook is inspired from the available online notebook [Ambient Seismic Noise Analysis](https://krischer.github.io/seismo_live_build/html/Ambient%20Seismic%20Noise/NoiseCorrelation.html) originally created by Celine Hadziioannou and Ashim Rijal. The original version of the notebook was modified by Léonard Seydoux (seydoux@ipgp.fr) in 2023 for the course "Scientific Computing for Geophysical Problems" at the [institut de physique du globe de Paris](http://www.ipgp.fr). It now includes a synthetic part to inspect the influence of the source distribution on the cross-correlation functions.

## Goals

This Jupyter notebook shows how to turn continuous records of ambient seismic noise into virtual seismograms using the theory of cross-correlation. The goal of this notebook is to reproduce part of the results obatined in the paper [_High-Resolution Surface-Wave Tomography from Ambient Seismic Noise_, by Shapiro, et al. (2005)](https://www.science.org/doi/10.1126/science.1108339).

###  Requirements

This notebook relies on Python libraries listed below. If you are running this notebook on the virtual machine provided for the course, you should have all the required libraries installed. If you are running this notebook on your own machine, you will need to install the following libraries:

- [ObsPy](https://github.com/obspy/obspy/wiki), an open-source project that provides a Python framework for processing seismological data. It provides parsers for standard file formats, clients to access data centers, and signal processing routines that allow the manipulation of seismological time series.
- [NumPy](https://numpy.org), an open-source project aiming to enable numerical computing with Python.
- [SciPy](https://numpy.org), fundamental algorithms for scientific computing in Python.
- [Matplotlib](https://matplotlib.org), a comprehensive library for creating static, animated, and interactive visualizations in Python.

If unavailable on your machine, you can use the [Anaconda](https://www.anaconda.com) package manager to install these libraries. Please ensure you can run the following cell without any error before proceeding further.

## Material

The notebook named [cross_correlation.ipynb](cross_correlation.ipynb)
presents the different concepts and examples with exercices to be completed. The solution to this notebook is provided in the notebook named [cross_correlation_solution.ipynb](cross_correlation_solution.ipynb). The data needed to run the notebook are downloaded directly in the notebook, which implies that an internet connection is required to run the notebook.

## Notebook contents

1. Introduction
    1. Goals of this notebook
    2. Requirements
    3. Theoretical background
2. Numerical experiments
    1. Experimental setup
    2. Generation of the synthetic
    3. Cross-correlation of the synthetic 
    4. Cross-correlation of the synthetic seismograms for different azimuths
    5. Mixture of all sources
    6. Getting rid of cross terms
3. Application to real data
    1. Data download
    2. Data inspection
    3. Data processing
    4. Cross-correlation of the data
    5. Comparison of the cross-correlation with an earthquake
4. Perspectives
    