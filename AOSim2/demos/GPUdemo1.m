Nloops = 5;
N = 250;
GFX = false;

gpu = gpuDevice(1);
gpu.reset;

PS = AOScreen(512);
PS.name = 'Phase Screen';
PS.describe;
CUBE = zeros([512 512 N]);

[X,Y] = PS.COORDS;
colormap(gray(256));

for USE_GPU=[false true]
    PS.grid(single(PS.grid_));
    PS.gpuify(USE_GPU);
    
    for m=1:Nloops
        tic;
        PS.make;
        fprintf('Loop %d: %s make took %.4f secs, ',m,PS.describe,toc);
        
        tic;
        s=1;
        for n=1:N
            s = s * .99;
            CUBE(:,:,n) = gather(PS.interpGrid(s*X,s*Y));
        end
        fprintf('Zooming took %.4f secs.\n',toc);
        
        if(GFX)
            animLoop(CUBE,1);
        end
    end
end

clear CUBE

