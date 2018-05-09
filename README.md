# wrfclass
Initial python class that will plot WRF files, both deterministic and ensemble output. Basic variables at the surface and aloft are available currently with more to be added as the class progresses. These functions are combinations of ones I have created myself and found from various sources (e.g. Luke Madaus). Currently, functions support the plotting of:
```
Surface
=======
precipitable water
composite and simulated reflectivity
dewpoint, temperature, theta e, theta v, virtual temp, relative humidity
wind and basic combo-surface plots
Aloft
=====
Vorticity
Geopotential height
Temperature, relative humidity
Other
=====
Brightness temperature
Melting Level
Shear of two levels

```

Functions to add:
```
thickness
accumulated precipitation
precip probability (ensemble)
CAPE/CINH
LCL
Bunkers motion
skewt/sounding capability
cross section capability 
precip type
Red Flag Index
```
Convective Plots:
```
ML CAPE/CINH
SB CAPE/CINH
0-1km SRH/0-6km SRH
0-1km EHI/0-3km EHI
updraft helicity
max 10-m ws
SCP
STP
```

Ensemble Plots:
```
Probabilities (UH and Reflectivity)
Paintball
Mean and Spread
```

Vertical interpolation in both pressure and height are available to produce plots of variables aloft and those that consist of two vertical levels (e.g. 0-6km shear).

