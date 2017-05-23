function [OTF,TT,CENTER0] = computeOTF(PSF,CENTER)
% function [OTF,TT(disabled),CENTER0] = computeOTF(PSF,CENTER)
 
% CENTER0 = round(size(PSF)/2);
% [X1,X2,R] = mkImageCoords(PSF,1,CENTER0);
% Rotf = 20;
% Rtest = 17;

if(nargin<2)
    CENTER = findPeakCCD(PSF);
    % CENTER = CENTER0;
end

% CHEBY = 100;
% win = chebwin(65,80); % vert
% win1 = chebwin(size(PSF,1),CHEBY); % vert
% win2 = chebwin(size(PSF,2),CHEBY); % vert
% WIN = win1*win2';
% WIN = 1;

% SELECT = R(:)<Rtest;
% UV = [X1(:) X2(:)];

OTF = circshift(PSF,1-CENTER);
% OTF = circshift(WIN.*PSF,1-CENTER);
OTF = fft2(double(OTF));
% OTF = circshift(OTF,CENTER0);
OTF = fftshift(OTF);

% [OTF,TT] = alignOTF(OTF,UV,SELECT);

TT = []; 
% fprintf('WARNING: TT is disabled.\n');

return;
end
