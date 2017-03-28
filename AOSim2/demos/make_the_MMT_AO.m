%% How to build and program the AOSim2 MMT AO model.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.
% 20090504 JLCodona: Had to add a Seg.make line.  This is a bug that needs
% to be found.  Sorry.

%% Start clean...
% close
% clear classes

%% Load in some definitions...

% This limits the reconstructor probe spatial frequencies during
% programming.  This is to compensate for my different mode ordering so
% that when we say "56 modes" we get them all, instead of tippish and
% tiltish being way up like modes 83 and 84 and therefore not getting
% included at all.  This is a tricky point and has to have been dealt with
% in the MMT reconstructor design.  All roads lead back to GUIDO BRUSA!

% 3  6 10 15 21 28 36 45 55 66 78 91 105 120
% MAX_MODES = 56;
MAX_MODES = 120;

D = 6.5;
secondary = 0.1 * D;
aa = 0.01;

SPACING = 0.01;

% PUPIL_DEFN = [
%    0 0 D         1 aa 0 0 0 0 0
%    0 0 secondary 0 aa 0 0 0 0 0
%    0 0 spider   -2 aa 4 0 0 0 0];

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa 0 0 0 0 0 ];

% load data/PMMT.mat

Seg = AOSegment;
Seg.name = 'MMT Primary';
Seg.pupils = PUPIL_DEFN;
Seg.spacing(SPACING);
Seg.make;

clf;
% Seg.touch.make.show;
A = AOAperture();
A.name = 'MMT';
A.spacing(SPACING);
A.addSegment(Seg);
A.show;
colormap(gray);

%% Grab the actuator coordinates from my hexapolar DM design.  There are
% pix in the data/pix directory.  I have higher and lower actuator density
% options already designed as well.
load data/MMT_DM336_Actuators.mat ACT BAD

%% Make a DM with an OPD grid matched to the Aperture A...
DM = AODM(A);
DM.name = 'MMT DM336';

%% Add in the actuators.
DM.addActs(ACT*6.5/2,1,A);

% Mark BAD actuators
%DM.disableActuators(BAD);
%DM.disableActuators(MMT_BADACTS_NO_CURRENT);
%DM.disableActuators(MMT_BADACTS_NO_POSITION);

% Specify the boundary conditions... 
% DM.defineBC(5,8); % A circle of 8 null points at 5m radius.
DM.defineBC(5,8); % A circle of 8 null points at 5m radius.
DM.plotRegions; daspect([1 1 1]); drawnow;

%% Build the Shack-Hartmann WFS.

% WFS = AOWFS(A,D/12);
WFS = AOShackHartmann(A);
% WFS.defineSubApsByAperture(A,D/12);
WFS.defineSubApsByAperture(A,D/12,[1 1]*D/12/2); % Use the minimum number of lenslets.
WFS.name = 'MMT Shack-Hartmann WFS';
A.show; 
WFS.plotSubAps('r');
% WFS.quiver(); 
drawnow; % Show them.

%% Now for some real work.  Building the RECONSTRUCTOR...
RECON = AOReconstructor(A,DM,WFS);

% Now program this crazy thing.
% We can look at the singular values and choose things manually later.  
% For now, we will make default assumptions.  You can rebuild it quickly
% later.  The MMT currently runs with 56 modes corrected.

OWD = sqrt(MAX_MODES/pi);
% RECON.program(D,6*sqrt(2)); % Use Fourier modes. OWD is ~6 lambda/D for programming.
RECON.zprogram(D,11);  % program using Zernikes.
% RECON.dhprogram(D,10); % program using disk harmonics.
semilogy(RECON.s/RECON.s(1));

F = AOField(A);
F.padBy(64);
F.lambda = RECON.lambda;
[x,y] = coords(F);
% Look at the reconstructor for good "gut" feelings...
lim0 = min(RECON.RECONSTRUCTOR(:));
lim1 = max(RECON.RECONSTRUCTOR(:));

% for n=1:size(RECON.RECONSTRUCTOR,2)
%     DM.setActs(5*RECON.RECONSTRUCTOR(:,n)).touch.render;
    %F.planewave*A*DM;
    %     imagesc(x,y,F.interferometer(1),[0 3]);
    %     sqar;
    %     axis xy;
%     imagesc(x,y,DM.grid .* A.grid,[lim0 lim1]);colorbar;
%     daspect([1 1 1]);
%     drawnow;
% end

RECON.rebuild(500);
RECON.Nmodes
DM.nActs

% Mark BAD actuators
% DM.disableActuators(BAD);
% DM.disableActuators(MMT_BADACTS_NO_CURRENT);
% DM.disableActuators(MMT_BADACTS_NO_POSITION);

save data/MMTAO_Model_working A DM WFS RECON D % Paranoid.

% Now run something like the canned_GMT script, but don't load anything in
% first.  
