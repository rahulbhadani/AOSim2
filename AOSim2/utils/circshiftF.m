function Z = circshiftF(Z,OFFSET)

% function Z = circshiftF(Z,OFFSET)
% 
% This is like circshift, but works with sub-pixel shifts.
% Beware of complex artifacts.
% 
% JLCodona, 20080730

N1 = size(Z,1);
N2 = size(Z,2);

x1 = mkXvec(N1,1/N1);
x2 = mkXvec(N2,1/N2);

[X1,X2] = meshgrid(x1,x2);

Z = ifft2(exp(-2*pi*1i*(OFFSET(2)*X1+OFFSET(1)*X2)).*fft2(Z));
