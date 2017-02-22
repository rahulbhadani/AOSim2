classdef AODivField < AOField
	% AODivField Class: Divergent/Convergent version of AOField.
    % Use property DF.zFocus to set divergence.  
    % Whether it is diverging or converging depends on DF.z and
    % DF.direction.
    % zRef is where the grid is defined (things like grid spacing).
    % The coords, kcoords, spacing, extent, etc., depend on distance.
    % zFocus is a coordinate singularity.  Therefore, z too close to zFocus
    % generates a warning.  I hope to fix this someday.
    % In the meantime, I will generate a warning if you get too close to
    % the focus.  I may just revert to the collimated beam mode... I
    % haven't decided yet.  (20170222)
	
	%% Properties
	% Static Constants
	properties(Constant=true, SetAccess = 'private')
    end
	
	% Public properties
	properties(Access='public')
    end
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        zFocus = inf;   % Collimated beam.
        zRef = 1;       % Where the grid is defined.
    end
    
	% Private
	properties(Access='private' )
	end
	
	%% Methods
	methods
		function obj = AODivField(ref)
			obj = obj@AOField(ref);
        end

        function DF = setFocus(DF,zFocus)
            % Set everything needed for the field's overall focus.
            
            DF.zFocus = zFocus;
        end

        function DF = setRefZ(DF,zRef)
            % Set the grid definition distance.
            
            DF.zRef = zRef;
        end
        
    end
    
    %% static methods
    methods(Static=true)
        
    end
end
