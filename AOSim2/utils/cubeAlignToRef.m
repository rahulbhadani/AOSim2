function CUBE = cubeAlignToRef(CUBE,TEMPLATE,CENTER)

% CUBE = cubeAlignToRef(CUBE,TEMPLATE,[CENTER])
% 
% Center the frames to the TEMPLATE
% 
% 20120505 jlcodona: Added to the NECO toolbox.

SZ = size(CUBE);
if(nargin<3)
    CENTER = round(SZ(1:2)/2);
end

for n=1:size(CUBE,3)
    %IMG = conv2(single(CUBE(:,:,n)),W,'same');
    IMG = single(CUBE(:,:,n));
    PEAK = findPeakCCD(IMG);
    CUBE(:,:,n) = circshift(CUBE(:,:,n),CENTER-PEAK);
end
