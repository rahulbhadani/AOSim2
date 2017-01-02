function CUBE = cubeDivImg(CUBE,DIV)

% CUBE1 = cubeDivImg(CUBE,IMG)
% 
% This is a basic datacube operation.
% Divide each frame by IMG.  
% 
% BUGS: I don't check nothin'. 
%
% 20120504 jlcodona: Added to the NECO toolbox.
% 20131123 jlc: fixed this to multiply doubles.

for n=1:size(CUBE,3)
    IMG = double(CUBE(:,:,n));
    CUBE(:,:,n) = IMG ./ double(DIV);
end
