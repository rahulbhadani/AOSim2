function C = secondMoment(M,CENTER)

% C = secondMoment(M,{CENTER=[0 0]}); 
% Second moment of M in 2-D. (pixels)
% Returns an array with three numbers: [C11, C12, C22]

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

C = [sum(sum(X1.^2.*M)),sum(sum(X1.*X2.*M)),sum(sum(X2.^2.*M))]/M1;
