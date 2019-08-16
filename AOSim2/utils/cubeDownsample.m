function DCUBE = cubeDownsample(CUBE,N)

% DCUBE = cubeDownsample(CUBE,N)

N1 = size(CUBE,1);
N2 = size(CUBE,2);
N3 = size(CUBE,3);

DCUBE = zeros([N1 N2 floor(N3/N)],'single');

for n=1:(N3/N)
    FRAMES = (n-1)*N + (1:N);
    FRAMES(FRAMES<1) = [];
    FRAMES(FRAMES>N3) = [];
    
    DCUBE(:,:,n) = mean(CUBE(:,:,FRAMES),3);
end


