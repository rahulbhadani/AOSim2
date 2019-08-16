function CUBE = cubeMinusImg(CUBE,dIMG)

% CUBE = cubeMinusImg(CUBE,IMG)
% 
% This is a basic datacube operation.
% 
% 20120504 jlcodona: Added to the NECO toolbox.

for n=1:size(CUBE,3)
    IMG = double(CUBE(:,:,n));
    CUBE(:,:,n) = IMG - dIMG;
end
