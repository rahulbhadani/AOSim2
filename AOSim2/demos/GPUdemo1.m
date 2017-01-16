PS = AOScreen(512);
PS.name = 'Phase Screen';
PS.describe;

[X,Y] = PS.COORDS;
colormap(gray(256));

for GFX=[false true]
    PS.gpuify(false);
    PS.describe
    
    for m=1:5
        tic;
        PS.make;
        fprintf('Loop %d: %s make took %.4f secs\n',m,PS.describe,toc);
        
        tic;
        s=1;
        for n=1:250
            s = s * .99;
            G = PS.interpGrid(s*X,s*Y);
            if(GFX)
                imagesc([X(1) X(end)]*s,[Y(1) Y(end)]*s,G);
                sqar; axis xy;
                title(PS.describe);
                drawnow;
            end
        end
        fprintf('Just Zooming took %.4f secs\n',toc);
    end
    
    PS.gpuify(true);
    PS.describe
    
    for m=1:5
        tic;
        PS.make;
        fprintf('Loop %d: %s make took %.4f secs\n',m,PS.describe,toc);
        
        tic;
        s=1;
        for n=1:250
            s = s * .99;
            G = PS.interpGrid(s*X,s*Y);
            if(GFX)
                imagesc([X(1) X(end)]*s,[Y(1) Y(end)]*s,G);
                sqar; axis xy;
                title(PS.describe);
                drawnow;
            end
        end
        fprintf('Zooming (and drawing) took %.4f secs\n',toc);
    end
end