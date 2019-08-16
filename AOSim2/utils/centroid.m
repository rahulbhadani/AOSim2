function C = centroid(M,CENTER)

% C = centroid(M,[CENTER=[0 0]]): First moment of M in 2-D.
% M can be an image cube.

M = squeeze(M);
if(nargin<2)
    CENTER = [0 0];
end

C = [];

% M = M - mean(mean(M));
M1 = sum(M(:));
% M2 = mean(mean(abs(M).^2));

x1 = (1:size(M,1))-CENTER(1);
x2 = (1:size(M,2))-CENTER(2);

[X1,X2] = meshgrid(x1,x2);

C = squeeze([sum(sum(X1.*M)),sum(sum(X2.*M))]/M1);
