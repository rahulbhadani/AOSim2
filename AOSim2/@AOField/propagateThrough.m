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

    F.propagate2(F.z-Zlist(ZnextIndex))*ATMO.layers{ZnextIndex}.screen; % Fix WIND!!!

    if(isa(ATMO,'AOAtmo2'))
        F*ATMO.layers{ZnextIndex}.shadow; % Fix WIND!!!
        if(~isempty(ATMO.shadowUpdate))
            if(ATMO.shadowUpdate(ATMO,ZnextIndex,F))
                error('Error updating AOAtmo2 layer shadow.');
            end
        end
    end
    
    F.z = F.z - eps;
    
    %F.show; drawnow;
end

if(F.z > FinalZ)
    F.propagate2(F.z-FinalZ);
    %F.show;
end
