function CUBE1 = cubeMulImg(CUBE,IMG)

% CUBE1 = cubeMulImg(CUBE,IMG)
% 
% This is a basic datacube operation.
% Multiply each frame by IMG.  
% 
% BUGS: I don't check nothin'. 
%
% 20120504 jlcodona: Added to the NECO toolbox.
% 20131123 jlc: fixed this to multiply doubles.

CUBE1 = zeros(size(CUBE));

for n=1:size(CUBE,3)
    CUBE1(:,:,n) = double(CUBE(:,:,n)) .* double(IMG);
end
