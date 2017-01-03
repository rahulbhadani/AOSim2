function CUBE = cubeNormalize(CUBE,POWER_REF)

% CUBE = cubeNormalize(CUBE,[POWER_REF])
% 
% Make all the images have the same integral: 
% the original median value.
% 
% 20120504 jlcodona: Added to the NECO toolbox.

if(nargin<2)
    POWER_REF = 1;
end

SUMS = squeeze(sum(sum(CUBE)));
% SUM_ = median(SUMS);

for n=1:size(CUBE,3)
    %CUBE(:,:,n) = CUBE(:,:,n) / (SUMS(n)/SUM_);
    CUBE(:,:,n) = CUBE(:,:,n) / (SUMS(n)/POWER_REF);
end
