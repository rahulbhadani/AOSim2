classdef AOCoronagraph < AOSegment
    % AOCoronagraph class.
    % Works basically like an AOAperture.
    % Multiplying by an AOCoronagraph carries you to the Lyot plane.
    % 20151224: JLCodona.  UA.SO.CAAO.AOSim2
    
    properties
        REMAPPED    = []; % For use with PIAA.  If null, no PIAA.
        MASK        = [];
        APODIZER    = [];
        FPM         = [];
        LYOT        = []; % If null, then use the AOCoronagraph grid.
        Fnumber     = 1;
        D           = [];
        FL          = [];  % Focal length
       
        CENTROID = [];
        
        verbose     = false; % print debugging info.

        % I may want to protect this later.
    
        % This propagates the pupil field from the masked pixels to the
        % FPM.  
        Fmatrix = [];
        invFmatrix = [];
        s = [];
        antiLyot_U = []; 
        
        % switches
        DO_APODIZER = false;
        
        % AOField caches
        
        Fpp = [];
        Fpp_ = [];
        
        Ffp = [];
        Ffp_ = [];
        
        Flyot = [];
    end
    
    methods
        
        % Constructor
        function CORO = AOCoronagraph(varargin)
            CORO = CORO@AOSegment(varargin);
            CORO.name = 'Coronagraph';
            
            CORO.APODIZER = AOSegment(CORO);
            CORO.APODIZER.name = 'Pupil Apodizer';
            CORO.APODIZER.constant(1);
        end
        
        function CORO = setCentroid(CORO)
            % CORO = CORO.setCentroid();
            % Set the CENTROID property for determining the optical axis.

            CORO.CENTROID = CORO.centroid;
        end
        
        function CORO = setupFPM(CORO,FocalLength,lambdaRef,pixel_ld,FPM_MAX_ld)
            % CORO = setupFPM(CORO,FocalLength,lambda,pixel_ld,FPM_MAX_ld)
            % FocalLength: The coronagraph focal length to the FPM.
            % Fnumber: The Fnumber coming into the focal plane.
            % lambdaRef: The reference wavelength for computing the plate
            % scale.
            % pixel_ld: How big is the pixel in lambda/D
            % FPM_MAX_ld: How many lambda/D is the radius of the FPM region.
           
            if(isempty(CORO.MASK))
                CORO.setupMASK(1e-5);
            end
            
            CORO.FL = FocalLength;
            CORO.D = CORO.estimateD;
            CORO.Fnumber = CORO.FL/CORO.D;
            
            LAMBDA_D = CORO.Fnumber * lambdaRef;
            PIXEL = LAMBDA_D * pixel_ld;
            
            N = 2*ceil(FPM_MAX_ld / pixel_ld)+1;
            CORO.FPM = AOScreen(N);
            CORO.FPM.spacing(PIXEL);
            
            CORO.FPM.zero;
            CORO.FPM.name = 'FPM';
        end
        
        function CORO = setupMASK(CORO,thresh)
            % CORO = setupMASK(CORO,thresh)
            % Use grid>thresh defines boolean mask.
            
            CORO.MASK = logical(CORO.grid >= thresh);
        end

        function CORO = initAPODIZER(CORO)
            CORO.APODIZER = AOSegment(CORO);
            CORO.constant(1);
        end
        
        function CORO = initF(CORO,lambda)
            fprintf('Initializing the Fourier operator...\n');
            
            if(nargin<2)
                lambda = CORO.lambdaRef;
            end
            
            [Xp,Yp] = CORO.COORDS;
            [Xf,Yf] = CORO.FPM.COORDS;
            
            % Note that this doesn't include the measure factor.
            CORO.Fmatrix = exp((1i*pi/lambda/CORO.FL)*(Xf(:)*Xp(CORO.MASK)' ...
                + Yf(:)*Yp(CORO.MASK)'))...
                / sqrt(6*pi*sum(CORO.MASK(:)));
                %;
            % / pi;
            %* sqrt(sum(CORO.MASK(:))*prod(CORO.FPM.size));            
        end
        
        function CORO = mkInvF(CORO,thresh)
            fprintf('Building the pseudo-inverse...\n');
            if(nargin<2)
                thresh = 1e-8;
            end
            
            [U,S,V] = svd(CORO.Fmatrix,'econ');
            CORO.s = diag(S);
            CORO.invFmatrix = pseudoInv_rebuild2(U,CORO.s,V,CORO.s(1)*thresh);
            CORO.antiLyot_U = U;
        end        
        
        function CORO = PPtoLP(CORO,Fpp)
            % WARNING: This breaks the usual design model for AOSim2.
            % CORO.PPtoLP(Fpp)
            % All the resulting fields are kept in the CORO object.
            
            CORO.Fpp = Fpp;
            
            CORO.Fpp_ = CORO.Fpp.copy; % For post=aperture and apodizer.
            
            CORO.Ffp = AOField(CORO.FPM); % Pre-FPM.
            CORO.Ffp.name = 'Focal Plane Field';
            CORO.Ffp.lambda = Fpp.lambda;
            
            CORO.Ffp_ = CORO.Ffp.copy; % For post-FPM.
            
            CORO.Fpp_ * CORO; % pass through pupil.

            if(~isempty(CORO.APODIZER))
                CORO.Fpp_ * CORO.APODIZER; % Apodize.
            end

            if(isempty(CORO.Fmatrix))
                CORO.initF(CORO.Fpp_.lambda).mkInvF(1e-8);
            end
            
            %% Compute the focal plane fields
            
            CORO.Ffp.grid(...
                reshape(CORO.Fmatrix * CORO.Fpp_.grid_(CORO.MASK(:)),...
                CORO.FPM.size));
            CORO.Ffp_.grid(CORO.Ffp.grid); 
            
            CORO.Ffp_ * CORO.FPM;
            %CORO.Ffp_.grid(CORO.Ffp_.grid()-CORO.Ffp.grid());
            
            %% Compute the field in the Lyot plane.
            
            CORO.Flyot = CORO.Fpp.copy;  % Start with the relayed field.
            CORO.Flyot.name = 'Lyot Plane Field';

            Q = zeros(CORO.Flyot.size);
            Q(CORO.MASK(:)) = CORO.invFmatrix * (CORO.Ffp_.grid_(:)-CORO.Ffp.grid_(:));
            CORO.Flyot.grid(Q);
            
            CORO.Flyot.grid(CORO.Flyot.grid + CORO.Fpp_.grid );
            % Use this to update the FPM....
            
            CORO.Ffp_.grid(...
                reshape(CORO.Fmatrix * CORO.Flyot.grid_(CORO.MASK(:)),...
                CORO.FPM.size));
        end
        
        function COEFS = antiLyotCoefs(CORO,F)
            % COEFS = antiLyotCoefs([F])
            % Compute antiLyot coefs for Ffp_ unless F is given.
            
            if(nargin<2)
                COEFS = CORO.antiLyot_U' * CORO.Ffp_.grid_(:);
            else
                COEFS = CORO.antiLyot_U' * F.grid_(:);
            end
        end
        
        function rmFfp = rm_antiLyotModes(CORO,REMOVE,F)
            % rmFfp = CORO.rm_antiLyotModes(REMOVE,[F]);
            % Remove the selected modes from Ffp_ unless F is given.
            % Does not change any internal states.
            % REMOVE may be a boolean mask or a list of mode numbers.

            if(nargin<3)
                rmFfp = CORO.Ffp_.copy;
                rmFfp.name = 'post-FPM field without selected anti-Lyot modes';
                
                
                % Project Ffp_ onto anti-Lyot modes and remove selected set.
                rmFfp.grid(reshape(CORO.Ffp_.grid_(:) ...
                    - CORO.antiLyot_U(:,REMOVE) * ...
                    (CORO.antiLyot_U(:,REMOVE)' * CORO.Ffp_.grid_(:)),...
                    rmFfp.size));
                
            else
                rmFfp = F.copy;
                rmFfp.name = [F.name ' without selected anti-Lyot modes'];
                
                
                % Project Ffp_ onto anti-Lyot modes and remove selected set.
                rmFfp.grid(reshape(F.grid_(:) ...
                    - CORO.antiLyot_U(:,REMOVE) * ...
                    (CORO.antiLyot_U(:,REMOVE)' * F.grid_(:)),...
                    rmFfp.size));
            end
        end
    end
end

