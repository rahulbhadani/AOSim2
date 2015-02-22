function CUBE = cubeConv2(CUBE,KERNEL)

% CUBE = cubeConv2(CUBE,KERNEL)

for n=1:size(CUBE,3)
    CUBE(:,:,n) = conv2(CUBE(:,:,n),KERNEL,'same');
end
