% Testing Adaptive Optics
clear all;clf;close all;

addpath('/home/ivory/VersionControl/AOSim2/AOSim2');
addpath('/home/ivory/VersionControl/AOSim2/AOSim2/utils');
colormap(hot);

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.
PHOTONS_PER_EXPOSURE = 1e4;


% First create the pupil
x0 = 0; %X-coordinates of the center of the pupil
y0 = 0; %Y-coordinates of the center of the pupil
TRANSMISSION_TYPE = 1; % Tranmission Type is 1 for Mirror
WIND_SHIFT = [1 5]; % 1m/s of wind speed in x-direction, 5m/s in y-direction
D = 0.5; % Aperture diameter
SPACING = D/256; %Setting the number of pixels across the segment/pupil
SMOOTHING = SPACING/4;
PUPIL = [x0, y0, D, TRANSMISSION_TYPE, SMOOTHING, 0, 0, 0, 0, 0];

A = AOSegment;
A.spacing(SPACING);
A.name = 'AO Test';
A.pupils = PUPIL;
A.make;
A.show;

% Create a field
F = AOField(A);
% F.resize(256);
F.FFTSize = 512; %How many points are required compute FFT
F.lambda = lambda;

Theta_d = asin((F.lambda/D))*206265;
% Diffraction Angle in arc seconds. Sin theta = lambda/d

[x, y] = F.coords;  
% Maps a real life coordinate system to the pixels  based on spacing properies
% Since we pass A, a segment object to the field, it takes the values of
% spacing from segment object A.

%First set a field to be a planewave. 
% Then lets incident on the aperture which a single segment A in our case.
% 
figure;
F.planewave*A;
F.show;
colormap(hot);

%Create Point Spread Function
FOV = 5; %In arc seconds
PLATE_SCALE = Theta_d/10;
[PSF0, thx, thy] = F.mkPSF(FOV, PLATE_SCALE);
%mkPSF calculates Point Spread Function by doing Fast Fourier Transform
PSFMax = max(PSF0(:)); % Get max of PSF for normalizing the PSF
PSF0  = PSF0/PSFMax; %Normalize the PSF so that the brightest value is 1
figure;
imagesc(thx, thy, log10(PSF0),[-4 0]);
daspect([1 1 1]);
axis xy;
colorbar;
colormap(hot);

%% Now create some turbulence
Turbulence = AOScreen(2048);
%2048x2048 Grid size for Phase Screen,
%Make the phase screen bigger to get a good low-frequency behavior
Turbulence.name = 'My turbulence layer';
Turbulence.spacing(0.02);
Turbulence.setR0(0.15); %Fried parameter at 500nm
Turbulence.make;
[x_t, y_t] = Turbulence.coords; % Coordinate points on turbulence.

figure;
z = 1e3;
for i = 1:200
    subplot(2,2,1);
    Turbulence.show;
    Turbulence.shiftPixels([4 8]); % simulate wind.
    colormap(hot);
    drawnow;
    
    
    %We pass the plane through pase screen.
    F.planewave*Turbulence;
    F.propagate(z)*A;
    
    [PSF0, thx, thy] = F.mkPSF(FOV, PLATE_SCALE);
    %mkPSF calculates Point Spread Functieon by doing Fast Fourier Transform
    PSFMax = max(PSF0(:)); % Get max of PSF for normalizing the PSF
    PSF0  = PSF0/PSFMax; %Normalize the PSF so that the brightest value is 1
      subplot(2,2,2);
    imagesc(thx, thy, log10(PSF0),[-4 0]);
    colorbar;
    colormap(hot);
    drawnow;
end



