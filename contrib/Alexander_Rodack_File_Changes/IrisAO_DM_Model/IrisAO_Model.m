%% Clear Workspace
clear all;
clc;
close all;

%% Some Starting Comments
% This code makes a model of an IrisAO PTT111 type mirror.  It does not
% include the limits on stroke or tip/tilt angle that a real mirror has.
% The size of the segments is accurate to within 3 microns of a real
% mirror, and the piston, tip, and tilt functionality is correct (to the
% best of my knowledge right now, after lots of debugging).
%
% The scripts that are used in creating the mirror model and using it are
% all in the folder you found this one.  The ones that are important to
% look at should you want to adapt this script for your own purposes are:
%******************************
% makeIrisAODM.m
% IrisAOComputeZernPositions.m
% IrisAOPTT.m
%******************************
% These functions are what are used to make, and then update the shape of
% the mirror (the middle one is not required to be used, but can be helpful
% if you want to apply zernike polynomials).
%
% Please also note that this REQUIRES AOSim2 to be on your current pathway.
% The model is completely built in, and depends on the functionality of
% AOSim2.
%
% Finally, special thanks to Dr. John Codona for giving me the script that
% constructs the mirror.  Although it was debugged by me to remove some
% offset problems, it is an algorithm I would likely never have come up with.
% The application of piston, tip, and tilt is my own algorithm.  If you
% find a bug, or that it is not working accurately, please let me know what
% the problem you encountered is and I will work to fix it.  I am also in
% the process of writing a function to do the rendering to make it easier
% to port this model to separate scripts.
%
% ATR (email: atrodack@email.arizona.edu)



%% Set Initial Parameters and Flags
%Parameters (currently set to Specs from UAWFC Testbed)
segpitch = 606e-6;
magnification = 1;
lambda = AOField.VBAND;
GAMMA = 1;

%Flags
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed

verbose = true; %turns on/off plotting intermediate steps when adding PTT

plotsegshape = true; %turns on a subplot that has mirror segment displacement
%It is probably best to leave plotsegshape set to true, but you have the
%option if you feel like blindly trusting the rendering (which to the best
%of my knowledge is fine to do because I have debugged it for hours now...)
%verbose must be true to see the result of this

DEBUG = false; %allows checking of Segment Rendering....probably should can it
% WARNING: THIS SHOULD LIKELY ALWAYS BE FALSE. IT PLOTS IMAGES OF THE
% DISPLACEMENT MATRICES AS THEY ARE STORED IN VARIOUS PLACES TO ENSURE THE
% CORRECT GRID WAS PLACED AS THE MIRROR PTT IS BEING RENDERED.

%% Make the Mirror
DM = makeIrisAODM(magnification,verbose_makeDM);


% Grab the shape of a segment for future use (doesn't matter which one, at
% this point they are all the same)
mask = DM.segList{1}.Segment.grid_;

clc; %clears the command window of the text generated from adding segments to the AOAperture object

%Save the Coordinate vectors of the DM
[x,y] = DM.coords;

%Compute the extent of the Segment's Grid to Construct Aperture
extentx = abs(min(x)) + abs(max(x));
extenty = abs(min(y)) + abs(max(y));
fprintf('X Width = %0.4f mm\n',extentx * 10^3);
fprintf('Y Width = %0.4f mm\n',extenty * 10^3);

%% Make an Aperture and a Field
A = AOSegment(DM);
D = (extentx);
dx = (magnification *segpitch) / 100;
PNECO = [0 0 D 1 1.5*dx 0 0 0 0 0];
A.pupils = PNECO;
A.make;

DM.trueUp;
A.centerOn(DM); %Ensures the aperture is aligned with the DM. Previous
%versions of the IrisAO Model were found to be broken with this step due to
%an offset in the x-positioning of the mirror segments. That has been fixed

F = AOField(A);
F.name = 'IrisAO DM';
F.FFTSize = [2048 2048];
F.lambda = lambda;

%% Get Something to Put on Mirror
%There are all options of what can be sent to the mirror, in the syntax an
%actual IrisAO Mirror would like them, with a few minor differences. The
%actaul mirror wants its inputs in microns and milliradians, but this
%function wants meters and radians, and this includes a flag that tells the
%command to either overwrite the current settings (false) or bump them
%(true). If you are not bumping, only one should be uncommented at a time.
%See IrisAOPTT.m for more help. If using Zernikes, see
%IrisAOComputeZernPositions.m for help.

% Flatten the Mirror (remember if flattened, a 1 Angstrom piston is
% applied. Take care to always overwrite the PTTpos matrix when going from
% flat to positioned to avoid including this piston down the line)
% DM = IrisAOPTT(DM,1:37,zeros(37,1)*10^-6,zeros(37,1)*-10^-3,zeros(37,1)*10^-3,false);

%Apply only a Tip
% DM = IrisAOPTT(DM,1:37,zeros(37,1)*10^-6,ones(37,1)*10^-3,zeros(37,1)*10^-3,false);

% Apply a Tip and a Tilt
% DM = IrisAOPTT(DM,1:37,zeros(37,1)*10^-6,ones(37,1)*10^-3,ones(37,1)*-10^-3,false);

% Set a Random Mirror Shape
% DM = IrisAOPTT(DM,1:37,randn(37,1)*10^-6,randn(37,1)*-10^-3,randn(37,1)*10^-3,true);

% Set a Zernike Polynomial
Zernike_Number = [5];
Zernike_Coefficient_waves = [2];
PTTpos = IrisAOComputeZernPositions( lambda*1e6, Zernike_Number, Zernike_Coefficient_waves );
DM = IrisAOPTT(DM,1:37,PTTpos(1:37,1)*10^-6,PTTpos(1:37,2)*10^-3,PTTpos(1:37,3)*10^-3,false);


%% Do Piston Tip Tilt Separately from AOSegment
%Segment widths in the Model
segwidthx = 6.179e-4 * magnification;
segwidthy = 6.97e-4 * magnification;

%Segment Coordinate System
[xseg,yseg] = DM.segList{1}.Segment.coords;

%Build an AOGrid Object to do the rendering
M = AOGrid(DM);
M.name = 'IrisAO Mirror Grid';
M.spacing(DM.spacing);
M.zero;

%Build a second one to store displacement map if set by user
if plotsegshape == true
    MM = AOGrid(DM);
    MM.name = 'IrisAO Mirror Displacements';
    MM.spacing(DM.spacing);
    MM.zero;
end

%Go segment by segment.....
for ii = 1:length(DM.segList)
    % Get angles and pistons
    tip = DM.segList{ii}.tiptilt(1);
    tilt = DM.segList{ii}.tiptilt(2);
    piston = DM.segList{ii}.piston;
    
    % if mirror is set to flat, add a little piston so it keeps the
    % segments (sort of a bug, but it has to be done)
    if piston == 0
        if tip == 0
            if tilt == 0
                piston = 1e-10; %1 Angstrom of piston
            end
        end
    end
    
    % compute edge displacements
    tipheightmax = (segwidthx * sin(tip))/2;
    tiltheightmax = (segwidthy * sin(tilt))/2;
    
    % make vectors for displacement, across only defined segment pixels
    tipvec = linspace(-tipheightmax,tipheightmax,floor(segwidthx/DM.dx));
    tiltvec = linspace(-tiltheightmax,tiltheightmax,floor(segwidthy/DM.dx));
    
    %add zeros to vectors to bring them up to the length of the coordinate system
    sizing = DM.segList{ii}.Segment.size;
    tipdif = (sizing(1))-length(tipvec);
    tiltdif = (sizing(2))-length(tiltvec);
    
    if isequal(round(tipdif/2),tipdif/2)
        zerotip1 = zeros(1,tipdif/2);
    else
        zerotip1 = zeros(1,floor(tipdif/2));
    end
    if isequal(round(tiltdif/2),tiltdif/2)
        zerotilt1 = zeros(1,tiltdif/2);
    else
        zerotilt1 = zeros(1,floor(tipdif/2));
    end
    
    tipvec = horzcat(horzcat(zerotip1,tipvec),zerotip1);
    tiltvec = horzcat(horzcat(zerotilt1,tiltvec),zerotilt1);
    
    % Compute Tip/Tilt Grids
    [TIP,Y] = meshgrid(tipvec);
    [X,TILT] = meshgrid(tiltvec);
    
    % Compute Piston
    OPL = 2*piston; %mirror OPL
    
    %Compute Wavenumber
    k = (2*pi)/lambda;
    
    %Compute Segment Shape
    shapegrid = (TIP+TILT+OPL) .* mask;
    
    % Apply Shape to Segment
    if DEBUG == true
        nugrid = shapegrid; %Look at the actual displacement grid of the segments
    else
        nugrid = exp((1i*k).*(TIP+TILT+OPL)).*mask; %compute the complex exponetial instead
    end
    grid(DM.segList{ii}.Segment,nugrid);
    
    
    %% The Nasty bit of Rendering the Mirror
    N = AOGrid(1);
    % Grab the offsets so the code knows where to put the segment
    DM.segList{ii}.Segment.Offset = DM.segList{ii}.Offset;
    % Ensure the spacing is the same
    N.spacing(DM.segList{ii}.Segment.spacing);
    % Apply the offset
    N.Offset = DM.segList{ii}.Segment.Offset;
    % Set the Segment
    N.grid(DM.segList{ii}.Segment.grid_);
    % Add the Segment to the Mirror
    M = M + N;
    
    if plotsegshape == true
        grid(DM.segList{ii}.Segment,shapegrid);
        NN = AOGrid(1);
        % Grab the offsets so the code knows where to put the segment
        DM.segList{ii}.Segment.Offset = DM.segList{ii}.Offset;
        % Ensure the spacing is the same
        NN.spacing(DM.segList{ii}.Segment.spacing);
        % Apply the offset
        NN.Offset = DM.segList{ii}.Segment.Offset;
        % Set the Segment
        NN.grid(DM.segList{ii}.Segment.grid_);
        % Add the Segment to the Mirror
        MM = MM + NN;
    end
    
    % A snippet that allows the visual checking of what exactly is being
    % set to all the .grid_ properties.  This is where a lot of the error
    % was happening in the previous model, other than the computation that
    % does the projected piston thing John talks about. If the first 3
    % plots are the same, and the most recent segment added to the 4th
    % frame is the same as those 3, the rendering is working correctly
    if DEBUG == true
        figure(1)
        subplot(2,2,1)
        imagesc(xseg,yseg,(nugrid));
        sqar;
        axis xy;
        %     colorbar;
        %     caxis([-4e-6,4e-6]);
        title(sprintf('nugrid for \nSegment # %d\n',ii));
        subplot(2,2,2)
        imagesc(xseg,yseg,(DM.segList{ii}.Segment.grid_));
        sqar;
        axis xy;
        %     colorbar;
        %     caxis([-4e-6,4e-6]);
        title(sprintf('Segment.grid_ for \nSegment # %d\n',ii));
        subplot(2,2,3);
        imagesc(xseg,yseg,(N.grid));
        sqar;
        axis xy;
        %     colorbar;
        %     caxis([-4e-6,4e-6]);
        title(sprintf('N.grid_ for \nSegment # %d\n',ii));
        subplot(2,2,4)
        imagesc(x,y,(M.grid));
        sqar;
        axis xy;
        title(sprintf('Adding \nSegment # %d\n',ii));
        colorbar;
        %     caxis([-4e-6,4e-6]);
        colormap(jet);
        drawnow;
        pause(0.5);
    else
        if verbose == true
            if plotsegshape == true
                % The rendering is only partially trusted by the user. Plot
                % the complex PTT, the mirror shape as segments are added
                % (in terms of the displacement), and the mirror as the
                % light would see it (in terms of the complex exponential)
                figure(1)
                subplot(1,3,1)
                plotCAmpl(DM.segList{ii}.Segment.grid_,GAMMA);
                sqar;
                axis xy;
                title(sprintf('Complex PTT for \nSegment #%d\n',ii));
                subplot(1,3,2)
                imagesc(x,y,MM.grid);
                axis xy;
                sqar;
                title(sprintf('Adding Segment #%d\nMirror Segment OPL\n',ii));
                subplot(1,3,3)
                plotCAmpl(M.grid,GAMMA);
                sqar;
                axis xy;
                title(sprintf('Mirror as Light\nSee`s it'));
                drawnow;
                pause(0.1);
            else
                % The rendering is being trusted by the user. Just plot the
                % complex PTT and the mirror as segments are added
                figure(1)
                subplot(1,2,1)
                plotCAmpl(DM.segList{ii}.Segment.grid_,GAMMA);
                sqar;
                axis xy;
                title(sprintf('PTT for \nSegment # %d\n',ii));
                subplot(1,2,2)
                plotCAmpl(M.grid,GAMMA);
                sqar;
                axis xy;
                title(sprintf('Adding Segment #%d\nto IrisAO Mirror\n',ii));
                drawnow;
                pause(0.1);
            end
        end
    end
    
end

% Apply the rendered mirror shape to the DM variable
DM.combined.grid(M.grid);

%% Put a Field on it
F.planewave * DM * A;
if verbose == true
    figure(2);
else
    if DEBUG == false
        figure(1);
    else
        figure(2)
    end
end

F.show
colormap(gray);
drawnow;

%% Plot the PSF
if DEBUG == false
    THld = lambda/D * 206265; % Lambda/D in arcsecs.
    FOV = 50*THld; % FOV for PSF computation
    PLATE_SCALE = THld/5; % Pixel Size for PSF computation
    [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
    PSFmax = max(PSF(:));
    
    if verbose == true
        figure(3)
        colormap(gray);
    else
        figure(2)
    end
    
    % imagesc(thx,thy,PSF/PSFmax);
    imagesc(thx,thy,log10(PSF/PSFmax),[-3,0]);
    colormap(gray);
    axis xy;
    sqar;
    title('PSF for Current Mirror Positioning');
end
