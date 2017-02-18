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
        function ATMO = addLayer(ATMO,screen,alt)
            % ATMO = addLayer(ATMO,screen,alt)
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
                ATMO.layers{n}.screen.Offset = ATMO.layers{n}.Wind*ATMO.time;
                ATMO.layers{n}.shadow.Offset = ATMO.layers{n}.Wind*ATMO.time;
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
        
        %% path integration functions
        function opl = OPL_(ATMO,X,Y,z,BEACON)
            % opl = OPL_(ATMO,X,Y,z,BEACON)
            opl = zeros(size(X));
            
            if(nargin<5)
                BEACON = ATMO.BEACON;
            end
            
            % fprintf('DEBUG: OPL_ layers: ');
            for n=1:ATMO.nLayers
                % fprintf('%d ',n);
                if(~ATMO.layers{n}.ignore)
                    zLayer = ATMO.layers{n}.screen.altitude;
                    W = ATMO.layers{n}.Wind;
                    t = ATMO.time;
                    if(zLayer<=BEACON(3))
                        [XLayer,YLayer] = scaleCone(ATMO,X,Y,z,zLayer);
                        if(ATMO.layers{n}.screen.touched)
                            ATMO.layers{n}.screen.make;
                        end
                        opl_ = interpGrid(ATMO.layers{n}.screen,XLayer,YLayer);
                        opl_(isnan(opl_)) = 0;
                        opl = opl + opl_;
                        opl_ = interpGrid(ATMO.layers{n}.shadow,XLayer,YLayer);
                        opl_(isnan(opl_)) = 0;
                        opl = opl + opl_;
                    end
                end
            end
            %fprintf('\n');
            
            if(ATMO.GEOMETRY)
                opl = opl + ATMO.geomDistances(X,Y,z);
            end
        end
    
        %         function A = gpuify(A)
        %             for n=1:A.nLayers
        %                 A.layers{n}.screen.gpuify;
        %                 A.layers{n}.shadow.gpuify;
        %             end
        %         end
        
        
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