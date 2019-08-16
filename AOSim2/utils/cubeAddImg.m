function CUBE = cubeAddImg(CUBE,IMG)

% CUBE = cubeAddImg(CUBE,IMG)
% 
% This is a basic datacube operation.
% 
% 20120504 jlcodona: Added to the NECO toolbox.

for n=1:size(CUBE,3)
    CUBE(:,:,n) = CUBE(:,:,n) + IMG;
end
