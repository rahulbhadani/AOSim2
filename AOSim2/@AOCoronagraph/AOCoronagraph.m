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
    
    end
    
    methods
        
        % Constructor
        function CORO = AOCoronagraph(varargin)
            CORO = CORO@AOSegment(varargin);
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
            [Xp,Yp] = CORO.COORDS;
            [Xf,Yf] = CORO.FPM.COORDS;
            
            % Note that this doesn't include the measure factor.
            CORO.Fmatrix = exp((1i*pi/lambda/CORO.FL)*(Xf(:)*Xp(CORO.MASK)' ...
                + Yf(:)*Yp(CORO.MASK)'))...
                ;
            % / pi;
            %* sqrt(sum(CORO.MASK(:))*prod(CORO.FPM.size));            
        end
        
        function [Flyot,Ffp] = PPtoLP(CORO,Fpp)
            % WARNING: This breaks the usual design model for AOSim2.
            % [Flyot,Ffp] = PPtoLP(CORO,Fpp)
            % Note that Flyot is evaluated BEFORE the Lyot stop.
            
            Ffp = AOField(CORO.FPM);
            Ffp.name = 'Focal Plane Field';
            Ffp.lambda = Fpp.lambda;
            
            Fpp * CORO * CORO.APODIZER; % pass through pupil and apodizer.
           
            Ffp.grid(reshape(CORO.Fmatrix * Fpp.grid_(CORO.MASK(:)),CORO.FPM.size));
            Ffp * (Fpp.lambda*CORO.FL); % Normalizing is confusing.
            
            COMPLEX_FPMASK = exp(1i*Ffp.k*CORO.FPM.grid_(:)) - 1;
            Flyot = Fpp.copy;
            Flyot.name = 'Lyot Plane Field';
            Q = zeros(Flyot.size);
            
            Q(CORO.MASK(:)) = CORO.Fmatrix' * ( COMPLEX_FPMASK .* Ffp.grid_(:));
            Flyot.grid(Flyot.grid + Q);
        end
    end
end

