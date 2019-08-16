function CUBE = cubeMulImgSelf(CUBE,IMG)

% CUBE = cubeMulImgSelf(CUBE,IMG)
% 
% In-place version...
% This is a basic datacube operation.
% Multiply each frame by IMG.  
% 
% BUGS: I don't check nothin'. 
%
% 20120504 jlcodona: Added to the NECO toolbox.
% 20131123 jlc: fixed this to multiply doubles.

for n=1:size(CUBE,3)
    CUBE(:,:,n) = CUBE(:,:,n) .* IMG;
end
