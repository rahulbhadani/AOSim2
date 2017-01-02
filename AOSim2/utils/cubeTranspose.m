function CUBET = cubeTranspose(CUBE)

% CUBET = cubeTranspose(CUBE)
% 
% Transpose each image in the cube;
% 
% 20120504 jlcodona: Added to the NECO toolbox.

SIZE = size(CUBE);
SIZET(1) = SIZE(2);
SIZET(2) = SIZE(1);
SIZET(3) = SIZE(3);
CUBET = zeros(SIZET);

for n=1:size(CUBE,3)
    CUBET(:,:,n) = CUBE(:,:,n)';
end
