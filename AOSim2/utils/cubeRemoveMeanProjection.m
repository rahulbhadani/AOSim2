function CUBE = cubeRemoveMeanProjection(CUBE)

% CUBE = cubeRemoveMeanProjection(CUBE)

Q = mean(CUBE,3);
Q = Q(:);
Q = Q/norm(Q);

for n=1:size(CUBE,3)
    IMG = CUBE(:,:,n);
    IMG(:) = IMG(:) - Q*(Q'*IMG(:));
    CUBE(:,:,n) = IMG;
end
