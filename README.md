# Purpose
Routines built utilizing astropy coordinate frames and conversions 
that put a stellar group into a frame of reference by subtracting
off LSR motion and the average motion of the group. Plotting routines 
included. New ideas for converting to a reference frame welcome.
Example in analysis.py

# Input data
Data to be read is an output file from Chen's SNN output and data from apogee
to be joined on `gaia_source_id`. Header keys for input data must include:

```
chen gaia data: ra, dec, pmra, pmdec, parallax

cam apogee data: radial_velocity, gaia_source_id, source_id 

```
