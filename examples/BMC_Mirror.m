% An example of how to make a continuous face-sheet Deformable Mirror
%
% 20150412 ATRodack


%% BMC DM
%Parameters
lambda = AOField.RBAND; % Red light.
SPACING = 1e-5; % fine spacing

%DM Specs
nActs = 1020; %32x32 minus 4 in the corners
Max_Stroke = 1.5e-6;
Pitch = 225e-6;
sidelength = 3.5e-3;

%Coordinate System
xmin = -sidelength;
xmax = sidelength;
BMC_x = (xmin:SPACING:xmax);
ymin = -sidelength;
ymax = sidelength;
BMC_y = (ymin:SPACING:ymax);
[BMC_X,BMC_Y] = meshgrid(BMC_x,BMC_y);

%Construct Pupil Shape
BMC_pupil = ones([721,721]);
BMC_pupil = padarray(BMC_pupil,[250,250],0,'both');

%Construct the Pupil
A_BMC = AOSegment;
A_BMC.spacing(SPACING);
A_BMC.name = 'BMC Pupil';
A_BMC.grid(BMC_pupil);

%Make the Pupil a DM
BMC_DM = AODM(A_BMC);
[X,Y] = BMC_DM.COORDS;
BMC_DM.show;
colormap(gray);
title('BMC Pupil');
drawnow;
pause(2);

%% Create Actuator Locations
%Make Actuator Coordiante Space
actuator_x = xmin:Pitch:xmax;
actuator_y = ymin:Pitch:ymax;
[ACTUATOR_X,ACTUATOR_Y] = meshgrid(actuator_x,actuator_y);
%Turn off Actuators at corners
ACTUATOR_X(1,1) = 0; ACTUATOR_X(32,32) = 0; ACTUATOR_X(32,1) = 0; ACTUATOR_X(1,32) = 0;
ACTUATOR_Y(1,1) = 0; ACTUATOR_Y(32,32) = 0; ACTUATOR_Y(32,1) = 0; ACTUATOR_Y(1,32) = 0;
ACTUATOR_X(ACTUATOR_X==0) = [];
ACTUATOR_Y(ACTUATOR_Y==0) = [];
%Write the locations so AOSim2 will understand them
BMC_ACTS = zeros(1020,2);
BMC_ACTS(:,1) = ACTUATOR_X(:);
BMC_ACTS(:,2) = ACTUATOR_Y(:);

% Add Actuators
BMC_DM.addActs(BMC_ACTS);
% Define Boundary Conditions
BMC_DM = BMC_DM.defineBC(4e-3,5,'square');
% Plot Actuator Locations
BMC_DM.plotActuators;
title('BMC Actuator Locations');
pause(2);
% Plot Actuator Influence Regions
BMC_DM.plotRegions;
title('BMC Actuator Influence Regions');
pause(2);

%% Do something with the DM
% Flatten the Mirror
BMC_DM.flatten;


% Piston the mirror up by a micron
BMC_DM.bumpActs(ones(1020,1)*10^-6);
mesh(X,Y,BMC_DM.grid);
daspect([1 1 1e-3])
title('Pistoned BMC Mirror');
drawnow;