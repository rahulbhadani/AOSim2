function F = propagateThrough(F,ATMO,FinalZ)
  
% F = propagateThrough(F,ATMO,[FinalZ=0])
% 
% Propagate the field through an AOAtmo or AOAtmo2 using wave propagation.
% Starts at the z position in the AOField and ends at FinalZ or 0.
% 
% JLCodona

if(nargin<3)
    FinalZ = 0;
end

Zlist = ATMO.listHeights;

while(F.z>FinalZ)
    ZnextIndex = find((F.z>Zlist) & (Zlist>=FinalZ),1);
    if(isempty(ZnextIndex))
        break;
    end
    if(ATMO.verbosity>0)
        fprintf('AOField.propagateThrough z is now %.0f\n',F.z);
    end

    F.propagate2(F.z-Zlist(ZnextIndex))*ATMO.layers{ZnextIndex}.screen; % Fix WIND!!!

    if(ATMO.verbosity>0)
        clf;
        subplot(2,2,1);
        D = 1; %mean(F.extent);
        [x,y] = ATMO.layers{ZnextIndex}.screen.coords;
        imagesc(x,y,ATMO.layers{ZnextIndex}.screen.grid,[-1 1]*2e-6);
        axis xy; sqar; colorbar;
        setFoV(D);
        title(sprintf('phase screen %d',ZnextIndex));
    end
    
    if(isa(ATMO,'AOAtmo2'))
        F*ATMO.layers{ZnextIndex}.shadow; % Fix WIND!!!
        if(~isempty(ATMO.shadowUpdate))
            if(ATMO.shadowUpdate(ATMO,ZnextIndex,F))
                error('Error updating AOAtmo2 layer shadow.');
            end
        end
        
        if(ATMO.verbosity>0)
            subplot(2,2,2);
            %ATMO.layers{ZnextIndex}.shadow.show; setFoV(D);
            [x,y] = ATMO.layers{ZnextIndex}.shadow.coords;
            imagesc(x,y,ATMO.layers{ZnextIndex}.shadow.grid,[0 1]*1e-6);
            axis xy; sqar; colorbar;
            title(sprintf('shadow screen %d',ZnextIndex));
        end
    end

    if(ATMO.verbosity>0)
        subplot(2,2,3);
        F.plotC(4);
        setFoV(D);
        title('Complex Field');
        
        subplot(2,2,4);
        F.plotDex([-3 0]+3);
        setFoV(D);
        title('log Irradiance');
        
        drawnow;
        
        savePNG(sprintf('%d_atmo2_%.0f.png',ZnextIndex,ATMO.time*1e3),150);
    end
    
    F.z = F.z - eps;
    
    %F.show; drawnow;
end

if(F.z > FinalZ)
    F.propagate2(F.z-FinalZ);
    %F.show;
end
