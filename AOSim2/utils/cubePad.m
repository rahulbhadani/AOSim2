function CUBE1 = cubePad(CUBE,PADDING)

% CUBE1 = cubePad(CUBE,PADDING)
% 
% Pad a cube by PADDING on both sides in directions 1 and 2.
% 
% 20150522 JLCodona

TEST = padarray(CUBE(:,:,1),[1 1]*PADDING,'both');
SZ = size(TEST);

PADVAL = 0;

CUBE1 = zeros([SZ size(CUBE,3)]);

for n=1:size(CUBE,3)
    CUBE1(:,:,n) = padarray(CUBE(:,:,n),[1 1]*PADDING,PADVAL,'both');
end

