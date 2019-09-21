% Testing Adaptive Optics
clc;clear all;
clf;close all;
% 
 addpath('../AOSim2');
 addpath('../AOSim2/utils');
colormap(hot);

% lambda = AOField.RBAND; % Red light.
lambda = 1550e-9; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.
PHOTONS_PER_EXPOSURE = 1e4;

M = 1024/4;
% First create the pupil
x0 = 0; %X-coordinates of the center of the pupil
y0 = 0; %Y-coordinates of the center of the pupil
TRANSMISSION_TYPE = 1; % Tranmission Type is 1 for Mirror
WIND_SHIFT = [1 5]; % 1m/s of wind speed in x-direction, 5m/s in y-direction
D = .5; % Aperture diameter
SPACING = D/M; %Setting the number of pixels across the segment/pupil
SMOOTHING = SPACING/4;
PUPIL = [x0, y0, D, TRANSMISSION_TYPE, SMOOTHING, 0, 0, 0, 0, 0];

A = AOSegment;
A.spacing(SPACING);
A.name = 'AO Segment (An Aperture)';
A.pupils = PUPIL;
A.make;
A.show;

%%
% Create a field

% When selecting the FFT Size, you want to make sure to have a larger FFT
% gridspace than the size of the field. I have defined the variable M as
% the number of sample points in the field and set the FFT Size to 2M. One
% of the things to do to make sure parameters are chosen correctly is to
% look at PSF and see if it looks like an Airy disk and that the first null
% is at 1.22*Theta_d

F = AOField(A);
% F.resize(256*4);
F.FFTSize = M*2; %How many points are required compute FFT
F.lambda = lambda;

Theta_d = asin((F.lambda/D))*206265;
% Diffraction Angle in arc seconds. Sin theta = lambda/d
Theta_d*1.22
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
title('Planewave through the aperture');

%%
%Create Point Spread Function
FOV = 5; %In arc seconds
% This will result in the PSF being interpolated to a grid that has 
% samples over the diffraction limit
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


%% Create the AO System
% This portion is pulled from the make_kuiper_AO example

WFS_FPS = 500.;
AO_STARTTIME = 0.002;
gain=1;

STROKE = 5.5e-6; %% What is this?

DM = AODM(A);   %Creating a deformable mirror
DM.name = 'BMC 140 element MEMS';

xx = (1:12)*D/12;
xx = xx - mean(xx);
[X,Y] = meshgrid(xx,xx);
ACTS = [X(:) Y(:)];

DM.addActs(ACTS,1,A); %Add Actuators
DM.defineBC(D,8); % A circle of 8 null points at D radius.
figure;
A.show;
colormap(gray);
DM.plotRegions; daspect([1 1 1]);
drawnow;
title('Plane Wave on a BMC Deformable Mirror');


%% Build the Shack-Hartmann WFS.
figure;
WFS = AOWFS(A,D/12);
WFS.name = 'Omega SHWFS';
A.show; 
WFS.quiver(1); 
title('Shack-Hartmann Wave front Sensor');
drawnow; % Show them.

%% Now for some real work.  Building the RECONSTRUCTOR...
RECON = AOReconstructor(A,DM,WFS);

% Now program this crazy thing.
% We can look at the singular values and choose things manually later.  
% For now, we will make default assumptions.  You can rebuild it quickly
% later.  The MMT currently runs with 56 modes corrected.

RECON.adhocProgram(D/12*3);

RECON.zprogram(D,12);  % program using Zernikes.
RECON.rebuild(56).show;
DM.setActs(0);
colorbar;
title('Reconstructor');


%% Now create some turbulence
Turbulence = AOScreen(2048);
%2048x2048 Grid size for Phase Screen,
%Make the phase screen bigger to get a good low-frequency behavior
Turbulence.name = 'A Turbulence layer';
Turbulence.spacing(0.02);
Turbulence.setR0(0.15); %Fried parameter at 500nm
Turbulence.make;
[x_t, y_t] = Turbulence.coords; % Coordinate points on turbulence.
Turbulence.show;
colormap(hot);
colorbar;


%%
figure;
z = 1e2;
for i = 1:100
    subplot(2,2,1);
    Turbulence.show;
    Turbulence.shiftPixels([1 2]); % simulate wind.
    colormap(hot);
%     drawnow;
    
    F.planewave*Turbulence;
    F.propagate(z)*A;
    
    [PSF0, thx, thy] = F.mkPSF(FOV, PLATE_SCALE);
    %mkPSF calculates Point Spread Function by doing Fast Fourier Transform
    PSFMax = max(PSF0(:)); % Get max of PSF for normalizing the PSF
    PSF0  = PSF0/PSFMax; %Normalize the PSF so that the brightest value is 1
      subplot(2,2,2);
    imagesc(thx, thy, log10(PSF0),[-4 0]);
    colorbar;
    colormap(hot);
    title(['Uncorrected at Iteration, i = ', num2str(i)]);
%     drawnow;
    
    
    % APply AO
    t = i/WFS_FPS;
	
    WFS.sense(F*DM);
    if(t>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
        DM.setActs(0);
        DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes); 
        DM.removeMean.clip(STROKE);
        RECON.RECONSTRUCTOR ;
    end
    
    
    subplot(2,2,3);
    F.show();
    WFS.quiver(1);
    F*DM;
    [PSF0, thx, thy] = F.mkPSF(FOV, PLATE_SCALE);
    PSFMax = max(PSF0(:)); % Get max of PSF for normalizing the PSF
    PSF0  = PSF0/PSFMax; %Normalize the PSF so that the brightest value is 1
      subplot(2,2,4);
    imagesc(thx, thy, log10(PSF0),[-4 0]);
    colorbar;
    colormap(hot);
     title(['Corrected at Iteration, i = ', num2str(i)]);
    drawnow;
end


%We pass the plane through pase screen.

