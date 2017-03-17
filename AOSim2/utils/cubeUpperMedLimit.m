function CUBE = cubeUpperMedLimit(CUBE,sz,bias)

% CUBE = cubeUpperMedLimit(CUBE,sz,[bias])

if(nargin<3)
    bias = 0;
end

for n=1:size(CUBE,3)
    IMG = CUBE(:,:,n);
    IMGmed = medfilt2(IMG,[1 1]*sz);
    MASK = (IMG>IMGmed+bias);
    IMG(MASK) = IMGmed(MASK)+bias;
    CUBE(:,:,n) = IMG;
end
