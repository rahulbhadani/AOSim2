function CUBE = cubeLPF(CUBE,LIMIT,ORDER)

% CUBE = cubeLPF(CUBE,LIMIT,ORDER)
% 
% This is a basic datacube operation.
% 
% 20140329 jlcodona: Added to the OSU toolbox.

% x1 = fftshift(mkXvec(size(CUBE,1),1));
% x2 = fftshift(mkXvec(size(CUBE,2),1));
x1 = (mkXvec(size(CUBE,1),1));
x2 = (mkXvec(size(CUBE,2),1));
[X1,X2] = meshgrid(x1,x2);
R = sqrt(X1.^2 + X2.^2);
if(nargin<3)
    ORDER = 2;
end
FILTER = exp(-(R/LIMIT).^ORDER);

for n=1:size(CUBE,3)
    IMG = double(CUBE(:,:,n));
    IMG = ifft2(FILTER.*fft2(IMG));
    CUBE(:,:,n) = IMG;
end
