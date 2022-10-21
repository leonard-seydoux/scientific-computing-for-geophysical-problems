# Location of earthquake from real waveforms

> Homework delivered on Oct. 26, 2022. This should be given back before the next session. You can either go with a Jupyter notebook, or a latex report. The results and anwers to my questions should be emailed to me at seydoux@ipgp.fr.

## Read data

Using your expertise, and the knowledge your acquired during the first and second notebook session, read the seismograms downloaded for you under the `data/` folder. 

These seismograms correspond to 15 earthquakes located around Greece, between Italy and Turkey. There were downloaded using the script `download_earthquakes.py` and the file `data/earthquakes_greece.csv` obtained from the [USGS Search Earthquake Catalog tool](https://earthquake.usgs.gov/earthquakes/search/).

## Pick _P_ and _S_ arrival times

Using your newly acquired expertise as seismologist, identify the _P_ and _S_ wave arrivals visually (or automatically) on the waveform of each stations for an earthquake of your choice. Then report the arrival times in a csv file. 

Note that the automatic picking of phases can be done with ObsPy, without a garanteed quality. If you are curious, you can have a look at the [ObsPy's Trigger/Picker Tutorial](https://docs.obspy.org/tutorial/code_snippets/trigger_tutorial.html).

## Infer the earthquake epicenter

Assuming a distance law for $\tau_{SP}$, use the infered $\tau_{SP}$ to estimate the location of the 15 earthquakes. By comparing your obtained location with the trusted locations indicated by USGS (check the catalog file), discuss your obtained locations depending on the used stations, and propagation assumptions.