function [n,m] = Noll_to_Zernike(j)
% [n,m] = Noll_to_Zernike(j)
% 
% Convert sequential Noll index to n,m Zernike indices.
%
% j is the linear Noll coordinate, 
% n is the radial Zernike index 
% m is the azimuthal Zernike index.
%
% Code based on: 
% https://github.com/tvwerkhoven/libtim-py/blob/master/libtim/zern.py

if(j<1)
    error('Noll indices start at 1.');
end

n = 0;
m = 0;
j1 = j-1;
while (j1 > n)
    n = n + 1;
    j1 = j1 - n;
    m = (-1)^j * (mod(n,2) + 2*floor((j1+mod(n+1,2))/2));
end

% fprintf('The Zernike mode is [%d,%d].\n',n,m);
