function CUBE = cubeCircShift(CUBE,OFFSET)

% CUBE = cubeCircShift(CUBE,OFFSET)

for n=1:size(CUBE,3)
    CUBE(:,:,n) = circshift(CUBE(:,:,n),OFFSET);
end

