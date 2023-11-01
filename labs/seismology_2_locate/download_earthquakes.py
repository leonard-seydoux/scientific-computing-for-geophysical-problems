"""Download earthquakes waveforms.

This script downloads earthquake waveforms from the CSV file obtained from the
USGS earthquake search facility. 

Author: Leonard Seydoux (seydoux@ipgp.fr)

Date: Oct. 2022
"""

import obspy
import pandas as pd
import tqdm

from obspy.clients.fdsn import Client

FILEPATH_EARTHQUAKES = "data/earthquakes_greece.csv"
FILEPATH_WAVEFORMS = "data/earthquake_{}.mseed"
DATA_PROVIDER = "IRIS"
SEISMIC_NETWORK = "HL"
SEISMIC_STATION = "*"
CHANNEL = "*"
LOCATION = "*"
TIME_BEFORE_EVENT_SEC = 10
TIME_AFTER_EVENT_SEC = 600

# Read catalog
catalog = pd.read_csv(FILEPATH_EARTHQUAKES, parse_dates=[0])

# Speak out
print(f"Found {catalog.shape[1]} earthquakes from the catalog.")

# Connect to client
client = Client(DATA_PROVIDER)

# Loop over every event
waitbar = tqdm.tqdm(enumerate(catalog.time), desc="Downloading")
for index, time in waitbar:
    # Turn into obspy date
    starttime = obspy.UTCDateTime(time)
    waitbar.set_description(str(time))

    # Collect stream
    stream = client.get_waveforms(
        network=SEISMIC_NETWORK,
        station=SEISMIC_STATION,
        channel=CHANNEL,
        location=LOCATION,
        starttime=starttime - TIME_BEFORE_EVENT_SEC,
        endtime=starttime + TIME_AFTER_EVENT_SEC,
    )

    # Save stream
    stream.write(FILEPATH_WAVEFORMS.format(index))
