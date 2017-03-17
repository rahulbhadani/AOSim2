% Program to make a PHASE SCREEN for use with other programs.

fprintf('Phase Screen Writer.\n\n');

fprintf('Enter phase screen specifications...\n');
dx = input('Enter pixel spacing (m): ');
N = input('Enter grid size in pixels: ');
Cn2 = input('Enter turbulence structure constant Cn2 (m^(-2/3) at default wavelength 556 nm): ');
THICKNESS = input('Enter turbulence layer thickness (m): ');
LAMBDA = input('Enter the wavelength for the phase screen (nm): ')*1e-9;

PS = AOScreen(N);
PS.name = 'Phase Screen';
% PS.double_precision = 'true'; % Uncomment this to render in double precision.

PS.spacing(dx);
PS.TURBULENCE_MODEL = AOScreen.VON_KARMAN;  % Change this to one of the supported models.
PS.lambdaRef = AOField.VBAND;  % Set this to the wavelength where you set the Cn2 value.
PS.setCn2(Cn2,THICKNESS);

r0 = PS.r0;
fprintf('This should give a Fried length (r0) of %.2g m at reference wavelength %.0f nm.\n',r0,PS.lambdaRef*1e9);
fprintf('At wavelength %.0f nm the Fried length should be %.2g m.\n',LAMBDA*1e9,scaleR0(r0,PS.lambdaRef,LAMBDA));

PS.make;
% PS.show;  % Uncomment this to make a plot of the generated screen (wavefront deflection);

PHASE = PS.grid * 2*pi/LAMBDA;

OUTPUT_FILE = input('Enter the output filename: ','s');

save(OUTPUT_FILE,'PHASE','LAMBDA','dx','Cn2','THICKNESS','LAMBDA');
