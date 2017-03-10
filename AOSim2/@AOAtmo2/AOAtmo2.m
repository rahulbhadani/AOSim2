classdef AOAtmo2 < AOAtmo
    % AOAtmo2: AOAtmo class with a shadow phase screen for other physical models..
    %
    % See AOAtmo.
    
    % properties
    properties(Access='public')
        shadowUpdate = []; % function ref: error = shadowUpdate(atmo,nscreen,field)
    end
    
    methods
        % Constructor
        function ATMO2 = AOAtmo2(varargin)
            % note that I don't really support all of the AOScreen args in
            % AOAtmo*.  I just need the size of the output grid, so I will
            % throw away most of the details.  This should stop vararg
            % nesting.
            ATMO2 = ATMO2@AOAtmo(varargin{1});
            
            % Put in some default physics values.  (Not real units.)
            ATMO2.physics.dt = 1e-3;
            ATMO2.physics.csound = 340;
            ATMO2.physics.gain = 1e-3;
        end
        
        % Operations
        function ATMO = addLayer(ATMO,screen,alt,MASK)
            % ATMO = addLayer(ATMO,screen,alt,[MASK])
            % Add another atmosphere layer at altitude alt.
            % Also add a shadow screen for other physics.
            % The MASK is optional, but if included will be multiplied into
            % the field during interactions.  
            % The MASK should be an AOGrid of some type, for example, an
            % AOSegment.
            
            n = length(ATMO.layers)+1;
            L = struct; % start building the layer.
            L.name = sprintf('Layer %d:%s',n,screen.name);
            if(nargin>2)
                screen.altitude = alt;
            end
            
            screen.lambdaRef = ATMO.lambdaRef; % If this causes inconvenience, please tell JLCodona.
            
            L.screen = screen;
            L.shadow = screen.copy;
            L.shadow.name = ['Shadow of ' L.screen.name];
            L.shadow.zero;
            %L.shadow.policy.motionMode = 'wind'; 
            L.shadow.policy.motionMode = 'fixed';

            if(nargin>3) % there is a mask
                L.mask = MASK.copy;
            else
                L.mask = [];
            end
            
            L.Wind = [0 0];
            
            L.ignore = false;
            
            ATMO.layers{end+1} = L;
            ATMO.touch;
            
            if(ATMO.verbosity>0)
                fprintf('AOAtmo now has %d layers.\n',length(ATMO.layers));
            end
        end
        
        function ATMO = clearShadows(ATMO)
            for n=1:ATMO.nLayers
                ATMO.layers{n}.shadow.zero;
            end
            ATMO.touched = true;
        end
        
%         function ATMO2 = tick(ATMO2)
%             % Bump the clock by physics.dt
%             ATMO2.setObsTime(ATMO2.time+ATMO2.physics.dt);
%         end
        
        function ATMO = setObsTime(ATMO,t) % set observation time.
            %ATMO = setObsTime(ATMO,t) % set observation time.
            ATMO.time = t;
            
            for n=1:ATMO.nLayers
                % Always move the screen with the wind.
                ATMO.layers{n}.screen.Offset = ATMO.layers{n}.Wind*ATMO.time;
                
                % Move the shadow according to the policy.
                if(strncmp(ATMO.layers{n}.shadow.policy,'wind',4))
                    ATMO.layers{n}.shadow.Offset = ATMO.layers{n}.Wind*ATMO.time;
                end
                
                % Move the mask according to the policy.
                if(~isempty(ATMO.layers{n}.mask) && strncmp(ATMO.layers{n}.mask.policy,'wind',4))
                    ATMO.layers{n}.shadow.Offset = ATMO.layers{n}.Wind*ATMO.time;
                end
            end
        end
        
       function ATMO = zero(ATMO)
            % ZERO: Set all grids of an AOAtmo2 to zero.
            
            for n=1:ATMO.nLayers
                ATMO.layers{n}.screen.zero;
                ATMO.layers{n}.shadow.zero;
            end
            
            ATMO.grid_(:) = 0;
            ATMO.fftgrid_ = [];
        end        
        
        function g = grid(ATMO,nugrid)
            % g = grid(ATMO,nugrid)
            if(ATMO.nLayers == 0)
                g = grid@AOScreen(ATMO,nugrid);
                return;
            end
            
            [X,Y] = COORDS(ATMO);
            g = ATMO.OPL(X,Y,ATMO.z);
            ATMO.grid_ = g;
        end

        function A = sumOPL(A)
            % A = A.sumOPL()
            % Sum all of the screens from the BEACON.
            
            [X,Y] = A.COORDS;
            A.grid(A.OPL_(X,Y,A.z,A.BEACON));
        end
        
        %% ray trace and physical optics integration functions
        function [opl,trans] = OPL_(ATMO,X,Y,z,BEACON)
            % [opl,trans] = ATMO.OPL_(X,Y,z,BEACON)
            
            opl = zeros(size(X));
            
            DO_TRANS = (nargout>1);
            if(DO_TRANS)
                trans = ones(size(X));
            end
            
            if(nargin<5)
                BEACON = ATMO.BEACON;
            end

            % fprintf('DEBUG: OPL_ layers: ');
            for n=1:ATMO.nLayers
                % fprintf('%d ',n);
                if(~ATMO.layers{n}.ignore)
                    zLayer = ATMO.layers{n}.screen.altitude;
                    if(zLayer<=BEACON(3))
                        [XLayer,YLayer] = scaleCone(ATMO,X,Y,z,zLayer);
                        if(ATMO.layers{n}.screen.touched)
                            ATMO.layers{n}.screen.make;
                        end

                        % Screen
                        if(ATMO.layers{n}.screen.periodic)
                            opl_ = ATMO.layers{n}.screen.interpGrid(XLayer,YLayer);
                        else
                            opl_ = ATMO.layers{n}.screen.interpGrid(XLayer,YLayer);
                        end
                        opl_(isnan(opl_)) = 0;                       
                        opl = opl + opl_;
                        
                        % Shadow
                        opl_ = ATMO.layers{n}.shadow.interpGrid(XLayer,YLayer);
                        
                        opl_(isnan(opl_)) = 0;
                        opl = opl + opl_;
                        
                        % Mask
                        if(DO_TRANS && ~isempty(ATMO.layers{n}.mask))
                            tx_ = ATMO.layers{n}.mask.interpGrid(XLayer,YLayer);
                            tx_(isnan(tx_)) = 0; % Opaque outside mask bounds.
                            trans = trans .* tx_;
                        end
                    end
                end
            end
            %fprintf('\n');
            
            if(ATMO.GEOMETRY)
                opl = opl + ATMO.geomDistances(X,Y,z);
            end
        end
        
        %%
        function ATMO = listLayers(ATMO)
            % ATMO.listLayers()
            % Print out info about the ATMO.
            
            for n=1:ATMO.nLayers
                fprintf('Layer %d: <%s> has a %dx%d phase screen with spacing %f m.\n',n,ATMO.layers{n}.screen.name,...
                    ATMO.layers{n}.screen.size,ATMO.layers{n}.screen.dx);
                fprintf('\theight = %f m\n',ATMO.layers{n}.screen.altitude);
                fprintf('\tthickness = %f m\n',ATMO.layers{n}.screen.thickness);
                fprintf('\tCn2 = %.2g (r0=%.3f for just this screen)\n',ATMO.layers{n}.screen.Cn2,ATMO.layers{n}.screen.r0);
                fprintf('\tWind = [%.1f %.1f] m/s\n',ATMO.layers{n}.Wind)
                fprintf('\tOffset=[%.3f %.3f] m\n',ATMO.layers{n}.screen.Offset);

                if(~isempty(ATMO.layers{n}.shadow))
                    fprintf('\tThere is a %dx%d shadow screen, with spacing %f\n',ATMO.layers{n}.shadow.size,ATMO.layers{n}.shadow.dx);
                end
                if(~isempty(ATMO.layers{n}.mask))
                    fprintf('\tThere is a %dx%d mask, with spacing %f\n',ATMO.layers{n}.mask.size,ATMO.layers{n}.mask.dx);
                end
                
                if(ATMO.layers{n}.screen.useGPU) 
                    fprintf('\tThe screen uses the GPU.\n');
                end
                if(ATMO.layers{n}.shadow.useGPU) 
                    fprintf('\tThe shadow screen uses the GPU.\n');
                end
                if(ATMO.layers{n}.mask.useGPU) 
                    fprintf('\tThe layer mask uses the GPU.\n');
                end
            end
            
            fprintf('\tThe total Fried Scale for a star would be %.3f m.\n',ATMO.totalFriedScale);
        end 
        
        function ATMO = showTransmission(ATMO)
            % Show the ray-traced transmission through the AOAtmo2.
            
            [x,y] = ATMO.coords;
            [X,Y] = ATMO.COORDS;
            
            [OPL,TX] = ATMO.OPL_(X,Y,0);

            imagesc(x,y,TX);sqar;
            axis xy;
            title(['Transmission of ' ATMO.name]);

            %drawnow;
        end
        
        %%
        function A = addGaussian(A,CENTER,amp,width,whichOnes)
            % AOAtmo2 = AOAtmo2.addGaussian(CENTER,amp,width,[whichOnes])
            % Add a Gaussian bump to each screen.
            % whichOnes is a list of bools saying whether to affect
            % [screen,shadow].
            % i.e. whichOnes = [true,false]; means add to screen but not to
            % shadow.
            % The default is whichOnes = [true,false].
            
            if(nargin<5)
                whichOnes = [true,false];
            end
            
            if(length(whichOnes) < 2)
                whichOnes(2) = ~whichOnes(1);
            end
            
            for n=1:A.nLayers
                if(whichOnes(1))
                   A.layers{n}.screen.addGaussian(CENTER,amp,width);
                end
                
                if(whichOnes(2))
                   A.layers{n}.shadow.addGaussian(CENTER,amp,width);
                end
            end
            %touch(A);
        end        

        %%
        function A = gpuify(A,useGPU)
            if(nargin<2)
                useGPU = true;
            end
            
            for n=1:A.nLayers
                A.layers{n}.screen.gpuify(useGPU);
                A.layers{n}.shadow.gpuify(useGPU);
            end
            
            if(useGPU)
                A.grid_ = gpuArray(A.grid_);
            else
                A.grid_ = gather(A.grid_);
            end
        end
    end
end
