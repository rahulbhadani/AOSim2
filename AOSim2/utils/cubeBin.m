function CUBE1 = cubeBin(CUBE,BINS)

% CUBE1 = cubeBin(CUBE,BINS)
% 
% Bin a cube
% 
% 20130911 JLCodona

if(numel(BINS) == 1)
    BINS(2) = BINS(1);
end

TEST = downsampleCCD(CUBE(:,:,1),BINS(1),BINS(2));
SZ = size(TEST);

CUBE1 = zeros([SZ size(CUBE,3)]);

for n=1:size(CUBE,3)
    CUBE1(:,:,n) = downsampleCCD(CUBE(:,:,n),BINS(1),BINS(2));
end

