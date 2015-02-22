% A demonstration of using AOSim2 to watch the evolution of various
% statistics beyond a Kolmogorov phase screen.
% Use the Kuiper 61" to view the field.
% 
% 20150221 JLCodona

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.

D = 1.54;
secondary = 14.5/100;

SPACING = 0.01;            % fine spacing makes nice pupil images but is really overkill.
aa = SPACING;              % for antialiasing.
% aa = 0.04;
spider = 0.0254;
% spider = 0.01;

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa/2 0 0 0 0 0
   0 0 spider   -2 aa 4 0 D/1.9 0 0
   ];

% Since this demo only uses one AOSegment, I will not use the AOAperture wrapper.  
% I will only use AOSegment.  This is fine for simple pupils.

A = AOSegment;
A.spacing(SPACING);
A.name = 'Kuiper 61inch Primary';
A.pupils = PUPIL_DEFN;
A.make;

clf;
colormap(gray);

A.show;
input 'Press ENTER to Continue...'

%% Make a Kolmogorov phase screen.
TURBULENCE = AOScreen(2048); % Make it big so we get good low-frequency behavior.
%TURBULENCE.lambdaRef = AOField.VBAND; %This is the default.

TURBULENCE.spacing(.02);
TURBULENCE.setR0(r0); 
TURBULENCE.make;

TURBULENCE.show;
input 'Press ENTER to Continue...'


%% Make an AOField object.

F = AOField(A);
F.resize(1024); % make it big to study the field before the pupil.
F.FFTSize = 1024; % Used to compute PSFs, etc.
F.lambda = lambda;

[x,y] = F.coords;

F.planewave*A;
F.show;

input 'Continue...'

% This adds a reference wave to the field and computes the intensity.
% imagesc(x,y,F.interferometer(1),[0 3]);
% sqar;
% axis xy;
% drawnow;
% input 'Continue...'

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

input 'Continue...'

PIXEL_RANGE = 513 + (-400:400);
xx = x(PIXEL_RANGE);
yy = y(PIXEL_RANGE);

PUPIL_PIXELS = abs(x)<0.6*D;

%% Now propagate from the screen to a large distance.  
% No aperture because we are studying the evolution of the field.

N1=2; N2=2;

% MAX_RANGE selected so that the Fresnel scale is 1/10 of the grid.
% GRID_SIZE = max(F.extent);
% ZMAX = (GRID_SIZE/50).^2/F.lambda;

% RANGES = [1:10,20:10:100,200:100:1000,2000:1000:1e4];
% RANGES = logspace(1,log10(ZMAX),100);
RANGES = logspace(1,6,200);
SELECT = 513+(-15:15); % points for the phasor plot.

for z=RANGES

    % Assume the phase screen is in the z= plane.
    
    % Start the field over each time rather than repeatedly propagating.  
    % It minimizes numerical issues..
    
    F.planewave*TURBULENCE;
    F.propagate(z);
    
    PSI = F.subGrid(PIXEL_RANGE,PIXEL_RANGE);
    IRRADIANCE = PSI .* conj(PSI); 
    
    subplot(N1,N2,1);
    imagesc(xx,yy,IRRADIANCE,[0 3]);
    daspect([1 1 1]);
    axis xy;
    colorbar;

    if(z<1000)
        title(sprintf('Irradiance z=%.1f m',z));
    else
        title(sprintf('Irradiance z=%.2f km',z/100));
    end
    subplot(N1,N2,2);

    %hold on;
    plot(F.subGrid(SELECT,SELECT),'k.','MarkerSize',1);
    %hold off;
    %drawCircles(1,[0 0],1,'r--');
    grid;
    xlim([-1 1]*3);
    ylim([-1 1]*3);
    daspect([1 1 1]);
    title('Complex Pupil Field');
    xlabel('real');
    ylabel('imag');
    
    F*A;
    subplot(N1,N2,3);
    % F.show;
    PSI = F.subGrid(PUPIL_PIXELS,PUPIL_PIXELS);
    plotComplex(PSI,3);
    axis xy;
    axis off;
    title('Complex Pupil Field');
    
    subplot(N1,N2,4);
    [PSF,thx,thy] = F.mkPSF(3,THld/5);
    imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
    daspect([1 1 1]);
    axis xy;
    colorbar;

    title('PSF');
    
    
    drawnow; 
end

