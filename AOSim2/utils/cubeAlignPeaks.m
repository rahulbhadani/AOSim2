function CUBE = cubeAlignPeaks(CUBE,CENTER)

% CUBE = cubeAlignPeaks(CUBE,[CENTER])
% 
% Center the frames' smoothed PEAKS on CENTER
% 
% 20120505 jlcodona: Added to the NECO toolbox.

SZ = size(CUBE);
if(nargin<2)
    CENTER = round(SZ(1:2)/2);
end
    
win = chebwin(25);
W = win*win';
W=W/sum(W(:));


for n=1:size(CUBE,3)
    IMG = conv2(single(CUBE(:,:,n)),W,'same');
    PEAK = findPeakCCD(IMG);
%     imagesc(CUBE(:,:,n));sqar;
%     hold on;
%     plot(PEAK(2),PEAK(1),'ro');
%     hold off;
%     drawnow;
    CUBE(:,:,n) = circshift(CUBE(:,:,n),CENTER-PEAK);
end
