function PS = make(PS)
% PS = make(PS)
% 
% New, extensible turbulent screen maker.  
% This method computes a new random realization of a screen.
% It uses the parameters in the AOScreen class to determine the model.
% 
% Supported turbulence models (specified by static class constants):
% 
% PS = AOScreen(N);
% PS.TURBULENCE_MODEL = model_number.
% 
% The model numbers are 
% AOScreen.KOLMOGOROV    -- pure power law turbulence spectrum.
% AOScreen.TATARSKI      -- KOLMOGOROV with an inner scale.
% AOScreen.VON_KARMAN    -- TATARSKI modified to include an outer scale (L0).
% AOScreen.MODIFIED_ATMO -- VON_KARMAN modified to include the Hill inner
%                           scale enhancement.
% All models now have a variable power law PHI_n index ALPHA.  Set this in
% the AOScreen before make()'ing it.
%
% I put the phase screen in GPU memory if grid_ is a gpuArray, but I always
% compute it in the CPU.  This is a fairly rare operation, so it is
% probably okay.  I'll review when I have more experience.

% Start with the envelope, the random part is common to all models.

% This effectively turns off the automatic machinery of turbulence screens.
% It actually just short-circuits the "touch" mechanism.  You can still
% "manually" use turbulence screens.  
% If you want to do that, set 
% PS.TURBULENCE_MODEL = AOScreen.DISABLED;
% during normal use, but before calling PS.make, do this...
%
% PS.TURBULENCE_MODEL = AOScreen.VON_KARMAN;
% PS.make();
% PS.TURBULENCE_MODEL = AOScreen.DISABLED;
%
% Now you will have a phase screen that doesn't respond to touches.  This
% is useful for shadow screen models.

if(PS.TURBULENCE_MODEL == AOScreen.DISABLED)
    PS.touched = false;
    return;
end

USE_GPU = PS.useGPU;

[KX,KY] = PS.KCOORDS;
KR2 = KX.^2 + KY.^2;
if(USE_GPU)
    KR2 = gpuArray(KR2);
end

kappa0 = 2*pi/PS.L0;
kappam = 5.92/PS.inner_scale;

switch PS.TURBULENCE_MODEL
    case AOScreen.KOLMOGOROV
        PSD = 0.033 * PS.Cn2 * KR2.^(-PS.ALPHA/2);
        
    case AOScreen.TATARSKI
        PSD = 0.033 * PS.Cn2 * KR2.^(-PS.ALPHA/2);
        PSD = PSD .* exp(-KR2./kappam^2); % Inner scale
        
    case AOScreen.VON_KARMAN
        PSD = 0.033 * PS.Cn2 .* (kappa0^2+KR2).^(-PS.ALPHA/2); % Outer scale
        PSD = PSD .* exp(-KR2./kappam^2); % Inner scale

    case AOScreen.MODIFIED_ATMO
        kappal = 3.3/PS.inner_scale;
        K = sqrt(KR2);
        PSD = 0.033 * PS.Cn2 .* (kappa0^2+KR2).^(-PS.ALPHA/2); % Outer scale
        PSD = PSD .* (1 + 1.802*(K/kappal) - 0.254*(K/kappal).^(7/6)); % Hill bump.
        PSD = PSD .* exp(-KR2./kappal^2); % Modified Inner scale
        
    otherwise % Use the Von Karman model as a default.
        PSD = 0.033 * PS.Cn2 .* (kappa0^2+KR2).^(-PS.ALPHA/2); % Outer scale
        PSD = PSD .* exp(-KR2./kappam^2); % Inner scale
end

PSD(PS.FAXIS_PIXEL(1),PS.FAXIS_PIXEL(2)) = 0;  % Kill DC for zero mean.
PSD = PSD * PS.thickness* 1.25; % The 1.25 is a fudge factor.  I'll have to figure that out later.  

PS.zero.addNoise(1);

PS.grid(PS.grid.*sqrt(PSD));

PS.grid(PS.fft);
PS.grid(PS.grid*(24.6 / sqrt(prod(PS.spacing)) / sqrt(prod(PS.size))  ));

if(PS.LOW_FREQ_FIX) % Patch up the real part using the imag part.
    
    % TBD
    fprintf('NOTICE: The Low Freq Patch is not included in this version YET.\n');
    
end

PS.real; % Throw away the imaginary part.  

PS.touched = false;

if(USE_GPU)
    PS.gpuify;
end

