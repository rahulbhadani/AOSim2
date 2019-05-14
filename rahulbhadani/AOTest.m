% Testing Adaptive Optics
addpath('/home/ivory/VersionControl/AOSim2/AOSim2');
addpath('/home/ivory/VersionControl/AOSim2/AOSim2/utils');

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.
PHOTONS_PER_EXPOSURE = 1e4;

WIND_SHIFT = [1 5]; % 1m/s of wind speed in x-direction, 5m/s in y-direction
D = 1.54; % Aperture diameter
SPACING = 0.01;