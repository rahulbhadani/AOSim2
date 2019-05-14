% Display the  phase screen

addpath('../AOSim2')
addpath('../AOSim2/utils')
% make a Kolmogorov phase screen
TURBULENCE = AOScreen(2048); % Make it big so we get good low-frequency behavior.
TURBULENCE.setR0(0.15); % Set fried parameter
TURBULENCE.make;
TURBULENCE.show;


%%
A = AOSegment;
A.spacing(SPACING);
A.name = 'Kuiper 61inch Primary';
A.pupils = PUPIL_DEFN;
A.make;

clf;
colormap(gray);

A.show;

%%
% Create a PlaneWave

LENS = AOScreen(A);
LENS.name = 'Lens';
LENS.zero.addZernike(2,0,-F.lambda/8,D);

F = AOField(A);
F.lambda = AOField.VBAND; %Wavelength of the field
F.planewave*A*LENS;
F.show;

[x,y] = F.coords;

%%
F.FFTSize = 1024; % Used to compute PSFs, etc.
THld = F.lambda/D * 206265; % Lambda/D in arcsecs.

% 
F.planewave*A; % Just go through the pupil.
[PSF,thx,thy] = F.mkPSF(5,THld/4);
PSFmax = max(PSF(:)); % Save for normalizing.

PSF = PSF/PSFmax; % make the brightest value =1.

imagesc(thx,thy,log10(PSF),[-4 0]); 
axis square;
axis xy;
colorbar;