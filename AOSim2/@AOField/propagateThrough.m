function F = propagateThrough(F,ATMO,FinalZ)
  
% F.propagateThrough(ATMO,[FinalZ=0])
% Propagate the field through an AOAtmo or AOAtmo2 using wave propagation.
% Starts at F.z and ends at FinalZ (defaults to 0).

nudge = 1e-6; % How far to push beyond a screen (just for ineq, could be eps).

if(nargin<3)
    FinalZ = 0;
end

Zlist = ATMO.listHeights;

if(sign(FinalZ - F.z) ~= F.direction)
    fprintf('WARNING: AOField "%s" is already past the destination plane.\n',F.name);
    return;
end

SCREENS = find(isBetween(Zlist,F.z,FinalZ));
while(~isempty(SCREENS))
    [~,POS] = min(abs(Zlist(SCREENS)-F.z));
    ZnextIndex = SCREENS(POS);
    Znext = Zlist(ZnextIndex);
    
    if(F.verbosity>0)
        fprintf('AOField "%s": Propagating to layer %d: %gm --> %gm (%gm)\n',...
            F.name,ZnextIndex,F.z,Znext,Znext-F.z);
    end
    F.propagate2(-(Znext-F.z));
    
    F*ATMO.layers{ZnextIndex}.screen; 
    if(isa(ATMO,'AOAtmo2'))
        F*ATMO.layers{ZnextIndex}.shadow; % Fix WIND!!!
        if(~isempty(ATMO.shadowUpdate))
            if(ATMO.shadowUpdate(ATMO,ZnextIndex,F))
                error('Error updating AOAtmo2 layer shadow.');
            end
        end
        
        if(ATMO.verbosity>1)
            subplot(2,2,2);
            %ATMO.layers{ZnextIndex}.shadow.show; setFoV(D);
            [x,y] = ATMO.layers{ZnextIndex}.shadow.coords;
            imagesc(x,y,ATMO.layers{ZnextIndex}.shadow.grid,[0 1]*1e-6);
            axis xy; sqar; colorbar;
            title(sprintf('shadow screen %d',ZnextIndex));
        end
    end
    
    F.z = F.z + F.direction * nudge;
    
    SCREENS = find(isBetween(Zlist,F.z,FinalZ));
end

if(F.z ~= FinalZ)
    if(F.verbosity>0)
        fprintf('AOField %s: Propagating from %gm to %gm (%gm)\n',...
            F.name,F.z,FinalZ,FinalZ-F.z);
    end
    F.propagate2(F.z-FinalZ);
end
