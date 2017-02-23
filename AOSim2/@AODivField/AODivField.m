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
        function DF = AODivField(ref)
            DF = DF@AOField(ref);
        end
        
        function DF = setFocus(DF,zFocus)
            % Set everything needed for the field's overall focus.
            
            DF.zFocus = zFocus;
        end
        
        function DF = setRefZ(DF,zRef)
            % Set the grid definition distance.
            
            DF.zRef = zRef;
        end
        
        function g = grid(obj,nugrid,mask) % TODO: make this fancier.
            %% function g = grid(obj,nugrid,[mask])
            % Dynamic version of grid.
            % Setting the grid assumes it is already demodulated.
            
            if(nargin==1) % Just read out the value.
                g = obj.grid_;
            else % Set the value to the input...
                if(nargin==2) % the CLASSIC behavior...
                    nugrid = squeeze(nugrid);
                    nugrid = nugrid(:,:,1);
                    nugrid = squeeze(nugrid);
                    %if(size(obj)==size(nugrid))
                    if(numel(obj)==numel(nugrid))  %  Also allow vector assignments
                        %obj.grid_(:) = nugrid(:); % Does not force shape.
                        obj.grid_ = reshape(nugrid(:),size(obj.grid_)); % Does not force shape.
                        %obj.touch;
                    else
                        obj.resize(size(nugrid));
                        obj.grid_ = nugrid;
                        %obj.touch;
                        %error('different sized grid assignment not supported (yet)');
                    end
                else % This when a mask is specified...
                    obj.grid_(mask(:)) = nugrid(:);
                end
                g = obj; % Note that I return the object if setting the grid.
            end
        end
        
        function scale = scalingFactor(DF,z)
            % scale = DF.scalingFactor(z)
            
            if(nargin<2)
                z = DF.z;
            end
            if(DF.zFocus ~= DF.zRef)
                if(DF.zFocus ~= 0)
                    scale = abs((z/DF.zFocus-1)/(DF.zRef/DF.zFocus-1));
                else
                    scale = abs(z/DF.zRef);
                end
                
            else
                error('ERROR: The Reference and Focus distances are the same. Singularity.');
                scale = 1;
                return;
            end
        end
        
        function s = spacing(DF,varargin)
            % s = spacing(obj,varargin)
            
            switch length(varargin)
                case 0
                case 1
                    arg = varargin{1};
                    if(isscalar(arg))
                        DF.spacing_ = [1 1]*arg;
                    else
                        DF.spacing_ = [arg(1) arg(2)];
                        DF.fftgrid_ = [];
                    end
                    
                case 2
                    DF.spacing_ = [varargin{1} varargin{2}];
                    DF.fftgrid_ = [];
                    
                otherwise
                    error('too many arguments');
            end
            
            s = DF.spacing_ * DF.scalingFactor(DF.z);
        end
        
        function D = extent(obj,sz)
            % D = obj.extent(new_extent)
            % No argument returns the spatial extent of the grid.
            % If the argument is a scalar, the object will be resized to
            % contain the new_extent.
            % If the argument is an AOGrid, the new extent is changed to
            % contain the object.
            
            if(nargin>1)
                if(isscalar(sz))
                    obj.resize(ceil(sz./obj.spacing));
                else
                    if(isa(sz,'AOGrid'))
                        obj.resize(ceil(sz.extent./obj.spacing));
                    end
                end
            end
            
            D = size(obj.grid_) .* obj.spacing;
        end

        function [X,Y] = COORDS(DF)
            [X,Y] = COORDS@AOField(DF);
            X = DF.scalingFactor * X;
            Y = DF.scalingFactor * Y;
        end
        
        function [X,Y] = coords(DF)
            [X,Y] = coords@AOField(DF);
            X = DF.scalingFactor * X;
            Y = DF.scalingFactor * Y;
        end
        
        %%
        function DF = show(DF,RANGE)
            % AODivField.show([RANGE]);
            % Make an image of the AOGrid.
            
            if(nargin<2)
                DF.plotC;
                title(DF.describe,'FontSize',14);
                daspect([1 1 1]);
            else
                DF.plotC(RANGE);
                title(DF.describe,'FontSize',14);
                daspect([1 1 1]);
            end
        end
    end
    
    %% static methods
    methods(Static=true)
        
    end
end
