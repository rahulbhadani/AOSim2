function CUBE = cubeULMedLimit(CUBE,sz,bias)

% CUBE = cubeUpperMedLimit(CUBE,sz,[bias])

if(nargin<3)
    % A bias of 0 is just like setting it equal to the median.
    bias = 50;
end

for n=1:size(CUBE,3)
    IMG = CUBE(:,:,n);
    IMGmed = medfilt2(IMG,[1 1]*sz);
    
    MASK = (IMG>IMGmed+bias);
    IMG(MASK) = IMGmed(MASK)+bias;
    
    MASK = (IMG<IMGmed-bias);
    IMG(MASK) = IMGmed(MASK)-bias;
    
    CUBE(:,:,n) = IMG;
end
