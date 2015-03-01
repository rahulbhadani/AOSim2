% I am assuming that the AOSim2 data directory is here.
% A symlink is good enough.

LAMBDA = AOField.JBAND;

D = 6.5;

PMMT = [    0            0          6.5            1         0.03            0            0            0            0            0
            0            0         0.64            0         0.02            0            0            0            0            0
            ];

Seg = AOSegment;
Seg.name = 'MMT Primary';
Seg.pupils = PMMT;
Seg.make;

clf;
% Seg.touch.make.show;
A = AOAperture;
A.name = 'MMT';
A.addSegment(Seg);
A.show;
colormap(gray);

% Make a test phase screen...
WFE = AOScreen(512);
WFE.name = sprintf('Mirror Residual WFE');
WFE.spacing(0.02);
WFE.setR0(D); % Diffraction-limited in VBAND.
WFE.make;

thld = LAMBDA/D*206265;
FOV = 25 * thld;
PLATE_SCALE = thld/3;

FINGER = AOSegment(A)
FINGER.name = 'The Finger!';

[X,Y] = FINGER.COORDS;

% Just put a mask into the grid directly...
% FINGER.grid(~(Y<-2.5 & abs(X-1)<0.25));
FINGER.grid(exp(1i*pi/2*double(Y<-2.5 & abs(X-1)<0.25)));

F = AOField(A);
F.lambda = LAMBDA;
F.FFTSize = 1024;

N1 = 3; 
N2 = 2;

subplot(N1,N2,1);
F.planewave*A*WFE;
F.show
title('Baseline pupil field');

% First part for the dOTF...
PSF0 = F.mkPSF(FOV,PLATE_SCALE);
PEAK = findPeakCCD(PSF0)
OTF0 = ifftshift(fft2(circshift(PSF0,1-PEAK)));

subplot(N1,N2,2);
F*FINGER;
F.show
title('Modified pupil field');

% Second part for the dOTF...
PSF1 = F.mkPSF(FOV,PLATE_SCALE);
OTF1 = ifftshift(fft2(circshift(PSF1,1-PEAK)));

dOTF = OTF0 - OTF1;

subplot(N1,N2,N2+1);
imagesc(log10(normalize([PSF0 PSF1])),[-4 0]);
daspect([1 1 1]);
axis xy;
axis off;
title('Baseline and mod PSFs');

subplot(N1,N2,N2+2);
imagesc(abs([OTF0 OTF1]));
daspect([1 1 1]);
axis xy;
axis off;
title('Baseline and mod MTFs');

subplot(N1,N2,2*N2+1);
imagesc(abs(dOTF));
daspect([1 1 1]);
axis xy;
axis off;
title('abs(dOTF)');

subplot(N1,N2,2*N2+2);
%imagesc(angle(dOTF));
phasor = exp(1i*angle(dOTF));
phasor(1,1) = 0; % so autoscaling doesn't mess up.
plotComplex(phasor,2);
daspect([1 1 1]);
axis xy;
axis off;
title('Phase of dOTF');

