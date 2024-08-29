# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 12:01:21 2022

@author: wcp_w2111150
"""

# NOTE YOU NEED GEOPANDAS

import datetime
import matplotlib.pyplot as plt
from cht_meteo.meteo import MeteoSource, MeteoGrid
from cht_meteo.meteo import filter_cyclones_TCvitals

# Define a track
source = MeteoSource("coamps", "coamps_tc_forecast", "forecast",
                     delay=6)

parameters = ["wind", "barometric_pressure"]
x_range = [-100.0, -80.0]
y_range = [10.0, 30.0]
ds = MeteoGrid(name="coamps", source=source, parameters=parameters,
               path="d:\data\COAMPS\hurricane_maria_coamps_new_t4",
               x_range=x_range,
               y_range=y_range)

t0 = datetime.datetime(2017,9,17,0,0,0)
t1 = datetime.datetime(2017,9,22,0,0,0)
ds.collect(time_range=[t0, t1])

# New determine best-track from these files
tracks = ds.find_cyclone_tracks(pcyc=99999, dt=6)

# Filter cyclone based on TCvitals
tracks = filter_cyclones_TCvitals(tracks)

# Write them out as geojson at location
for k, track in enumerate(tracks):
    
    # Check the size - only if there are 3 timestamps
    if track.track.geometry.size > 3:

        # Save as cyc
        output = r'd:\data\COAMPS\maria_tracks' + str(k) + '.cyc'
        track.write_track(output, 'ddb_cyc')

        # Save as shapefile
        output = r'd:\data\COAMPS\maria_tracks' + str(k) + '.shp'
        track.track.to_file(output)


# Use the first track to make ensembles
from cht_cyclones.tropical_cyclone import TropicalCyclone, TropicalCycloneEnsemble
from datetime import datetime, timedelta
import fiona
tc = tracks[0]
tc.account_for_forward_speed()
tc.estimate_missing_values()
tc.include_rainfall = True

tc2         = TropicalCycloneEnsemble(name="test", TropicalCyclone=tc)
tc2.tstart  = datetime(2017,9,17,0,0,0)
tc2.tend    = datetime(2017,9,22,0,0,0)
tc2.compute_ensemble(100)

fname = r'd:\data\COAMPS\ensembles'
tc2.to_shapefile(fname)
#tc2.to_spiderweb(fname)
#tc2.make_figures(fname)

# Done
print('done!')