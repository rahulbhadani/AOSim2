classdef AOCoronagraph < AOSegment
    % AOCoronagraph class.
    % Works basically like an AOAperture.
    % Multiplying by an AOCoronagraph carries you to the Lyot plane.
    % 20151224: JLCodona.  UA.SO.CAAO.AOSim2
    
    properties
        REMAPPED    = []; % For use with PIAA.  If null, no PIAA.
        PPMASK        = []; % Pupil Mask for reducing matrix rank.
        APODIZER    = [];
        FPM         = []; % Focal plane Phase mask
        FPTM        = []; % Focal plane Transmission mask
        
        LYOT        = []; % If null, then use the AOCoronagraph grid.
        Fnumber     = 1;
        D           = [];
        FL          = [];  % Focal length
       
        CENTROID = [];
        
        verbose     = false; % print debugging info.

        % I may want to protect this later.
    
        % This propagates the pupil field from the masked pixels to the
        % FPM.  
        
        % The Fourier matrix that carries Fpp_ to Ffp.  
        % I plan on having it use the full extent of Fpp_ and Ffp.  
        % Use the masks to reduce the domain and range.  
        % The pupil mask is imprinted on Fpp_ so it doesn't have to be in
        % the Fmatrix.
        Fmatrix = []; 
        
        % I don't think I need this anymore.
        invFmatrix = [];
        s = [];
        antiLyot_U = []; 
        
        % switches
        DO_APODIZER = false;
        
        % AOField caches
        
        Fpp = [];  % reference to the input pupil plane field.
        Fpp_ = []; % This is the field after apodization. PIAA remaps the pupil so this may be morphed.
        
        Ffp = [];  % Reference to the field in the focal plane, before the FPM.
        Ffp_ = []; % Field after the FPM.  Now probably useless.

        Flyot = []; % The field at the Lyot plane, but before any Lyot Stop is applied.
        
        Fscience = []; % The field at the final focal plane.

        % SELECTION MASKS
        
        FPM_ASSIGNED = [];
                
    end
    
    methods
        
        % Constructor
        function CORO = AOCoronagraph(varargin)
            CORO = CORO@AOSegment(varargin);
            CORO.name = 'Coronagraph';
            
            CORO.initAPODIZER();
            CORO.setupPPMASK(1e-5);
        end
        
        
        %% 
        
        function CORO = setupPPMASK(CORO,thresh)
            % CORO = setupPPMASK(CORO,thresh)
            % Use grid>thresh defines PUPIL boolean mask.
            
            CORO.PPMASK = logical(CORO.grid >= thresh);
        end

        function CORO = initAPODIZER(CORO)
            % CORO.initAPODIZER();
            CORO.APODIZER = AOSegment(CORO);
            CORO.APODIZER.name = 'Pupil Apodizer';
            CORO.APODIZER.constant(1);
        end
        
        function CORO = setupFPM(CORO,FocalLength,lambdaRef,pixel_ld,FPM_MAX_ld)
            % CORO = setupFPM(CORO,FocalLength,lambda,pixel_ld,FPM_MAX_ld)
            % FocalLength: The coronagraph focal length to the FPM.
            % Fnumber: The Fnumber coming into the focal plane.
            % lambdaRef: The reference wavelength for computing the plate
            % scale.
            % pixel_ld: How big is the pixel in lambda/D
            % FPM_MAX_ld: How many lambda/D is the radius of the FPM region.
           
            % Derive and estimate params
            CORO.lambdaRef = lambdaRef;
            CORO.FL = FocalLength;
            CORO.D = CORO.estimateD;
            CORO.Fnumber = CORO.FL/CORO.D;
            
            LAMBDA_D = CORO.Fnumber * lambdaRef;
            PIXEL = LAMBDA_D * pixel_ld;  % Fp plate scale.
            
            %% Make Focal Plane Masks
            N = 2*ceil(FPM_MAX_ld / pixel_ld)+1;
            CORO.FPM = AOScreen(N);
            CORO.FPM.spacing(PIXEL);
            
            CORO.FPM.zero;
            CORO.FPM.name = 'FPM OPD';
            
            CORO.FPTM = AOSegment(CORO.FPM);
            
            CORO.FPTM.constant(1.0);  % All pixels transmit.
            CORO.FPTM.name = 'FPM Transmission';
            
            CORO.FPM_ASSIGNED = false(CORO.FPM.size); % Nothing is assigned.
            
        end
        
        %%
        function CORO = initF(CORO,lambda)
        % CORO.initF(lambda);
        % Compute the operator that takes Fp_ from PP to FP before FPM.
        % NOTE: I use the FULL extent of Fp_ and Ffp for Fmatrix.
        % You may want to downselect rows and columns later with boolean
        % filters.
        % TODO: Fix the normalization so this is correct in Fourier Optics.
        
            fprintf('Initializing the Fourier operator...\n');
            
            if(nargin<2)
                lambda = CORO.lambdaRef;
            end
            
            [Xp,Yp] = CORO.COORDS;
            [Xf,Yf] = CORO.FPM.COORDS;
            
            % Note that this isn't the proper measure factor.
            CORO.Fmatrix = exp((1i*pi/lambda/CORO.FL)*...
                ( Xf(:)*Xp(:)' + Yf(:)*Yp(:)') );
                %* prod(CORO.spacing);
        end
        
        function CORO = mkInvF(CORO,thresh)
            % CORO = mkInvF(CORO,thresh)
            % I don't think I need this anymore.
            
            fprintf('Building the pseudo-inverse...\n');
            if(nargin<2)
                thresh = 1e-8;
            end
            
            [U,S,V] = svd(CORO.Fmatrix,'econ');
            CORO.s = diag(S);
            CORO.invFmatrix = pseudoInv_rebuild2(U,CORO.s,V,CORO.s(1)*thresh);
            CORO.antiLyot_U = U;
        end        
        
        %% Operators and tasks
        
        function PPtoFP(CORO)
            % CORO.PPtoFP()
            % Use Fmatrix to find Ffp starting from Fpp_.
            % Make sure Ffp_ is defined and properly initialized before
            % calling this.
            CORO.Fpp_ = CORO.Fpp.copy * CORO * CORO.APODIZER;
            CORO.Ffp.grid(CORO.Fmatrix * CORO.Fpp_.grid_(:));
        end
        
        function PPtoFPmask(CORO,usePPMask)
            % CORO.PPtoFP([usePPMask])
            % Use Fmatrix to find Ffp starting from Fpp_.
            % Make sure Ffp_ is defined and properly initialized before
            % calling this.
            % The optional flag usePPMask (defaults true) determines use of
            % PPMASK.
            
            if(nargin<2)
                usePPMask = true;
            end

            if(usePPMask)
                CORO.Ffp.grid(...
                    reshape(CORO.Fmatrix(:,CORO.PPMASK(:)) * CORO.Fpp_.grid_(CORO.PPMASK(:)),...
                    CORO.FPM.size));
            else
                CORO.Ffp.grid(reshape(CORO.Fmatrix * CORO.Fpp_.grid_(:),CORO.FPM.size));
            end
        end
        
        
        function CORO = PPtoLP(CORO,Fpp)
            % CORO.PPtoLP(Fpp)
            % 
            % This sets up the internal fields for whatever use.  
            % 
            % WARNING: This breaks the usual design model for AOSim2.
            % All the resulting fields are kept in the CORO object.
            
            CORO.Fpp = Fpp;
            
            CORO.Fpp_ = CORO.Fpp.copy; % For post=aperture and apodizer.
            CORO.Fpp_.name = 'post-APOD Field';
            
            CORO.Ffp = AOField(CORO.FPM); % Pre-FPM.
            CORO.Ffp.name = 'Focal Plane Field';
            CORO.Ffp.lambda = Fpp.lambda;

            CORO.Ffp_ = CORO.Ffp.copy; % For post-FPM.
            
            CORO.Fpp_ * CORO; % pass through pupil.

            if(~isempty(CORO.APODIZER))
                CORO.Fpp_ * CORO.APODIZER; % Apodize.
            end

            if(isempty(CORO.Fmatrix))
                CORO.initF(CORO.Fpp_.lambda);
                %CORO.initF(CORO.Fpp_.lambda).mkInvF(1e-8);
            end
            
            %% Compute the focal plane fields
            
            CORO.PPtoFP;
            
            CORO.Ffp_.grid(CORO.Ffp.grid); % This keeps the Ffp from being modified.
            CORO.Ffp_ * CORO.FPM * CORO.FPTM;
            
            %% Compute the field in the Lyot plane.
            
            CORO.Flyot = CORO.Fpp.copy;  % Start with the relayed field.
        end
        
        function CORO = FPtoLP(CORO,FPSELECT)
            % CORO.FPtoLP([FPSELECT])
            % Carry selected FP pixels to LP.
            % This currently goes to all points in LP, metrics should use PPMASK.
            % FPSELECT overrides FPM_ASSIGNED which is the default.
            
            if(nargin<2)
                FPSELECT = CORO.FPM_ASSIGNED;
            end
            
            psiF = CORO.Ffp.grid_(FPSELECT(:));
            psiF = psiF .* CORO.FPTM.grid_(FPSELECT(:));
            psiF = psiF .* exp(1i*CORO.Ffp.k * CORO.FPM.grid_(FPSELECT(:)));
            g = nan(CORO.Ffp_.size);
            g(FPSELECT(:)) = psiF;
            CORO.Ffp_.grid(g);
            
            psiL = CORO.Fmatrix(FPSELECT(:),:)' * psiF;
            %CORO.Flyot.grid(reshape(psiL,CORO.Flyot.size));
            CORO.Flyot.grid(psiL);
        end

        
        %% Coronagraph operations
        
        function CORO = response(CORO,ANGLE0)
            % CORO.response(ANGLE0)
            % Amplitude 1 response through the corograph.
            
            CORO.Fpp.planewave(1,ANGLE0);
            CORO.PPtoFP;

            CORO.Ffp_ = CORO.Ffp.copy;
            CORO.Ffp_ * CORO.FPTM;
            
            CORO.Flyot.grid(CORO.Fmatrix' * CORO.Ffp_.grid_(:));
            
            if(isempty(CORO.LYOT))
                CORO.Flyot * CORO;
            else
                CORO.Flyot * CORO.LYOT;
            end
            
            CORO.Fscience = CORO.Ffp.copy;
            CORO.Fscience.name = 'Science Cam';
            CORO.Fscience.zero;
            psi = CORO.Fscience.grid;
            %psi(CORO.PPMASK(:)) = CORO.Fmatrix * CORO.Flyot.grid_(:);
            psi(:) = CORO.Fmatrix * CORO.Flyot.grid_(:);
            CORO.Fscience.grid(psi);
        end
        
        
        %%
        function LPvectors = LyotContribs(CORO,FPSELECT,JustPP)
         % LPvectors = CORO.LyotContribs([FPSELECT],[JustPP])
         % FPSELECT defaults to ALL PIXELS in FPM
         % JustPP defaults to false.
         % The LP field vectors (no masks) for each pixel in FPSELECT.
         %
         % NOTE: This does not include the effects of the FPM.
         
            if(nargin<2)
                FPSELECT = true(CORO.FPM.size);
            end
         
            if(nargin<3)
                JustPP = false;
            end
            
            psiF = CORO.Ffp.grid_(FPSELECT(:));
            % psiF = psiF .* CORO.FPTM.grid_(FPSELECT(:));
            psiF = psiF .* exp(1i*CORO.Ffp.k * CORO.FPM.grid_(FPSELECT(:)));

            if(JustPP)
                LPvectors = CORO.Fmatrix(FPSELECT(:),CORO.PPMASK(:))';
            else
                LPvectors = CORO.Fmatrix(FPSELECT(:),:)';
            end
                
            for n=1:size(LPvectors,2)
                LPvectors(:,n) = LPvectors(:,n) * psiF(n);
            end
         
        end
        
        %% MAGIC!
        
        function CORO = reset(CORO,FPSELECT)
            % CORO.reset([FPSELECT]);
            % Reset the FPM to the starting value.
            
            if(nargin>1)
                if(size(FPSELECT) == CORO.FPM.size)
                    CORO.FPM_ASSIGNED = FPSELECT;
                else
                    fprintf('ERROR: FPSELECT must be the same size as the FPM');
                    return;
                end
            end
            
            CORO.FPM.zero;
        end
        
        %% First Incremental experiment.
        function CORO = algo1(CORO,FPSELECT)
            % CORO.algo1([FPSELECT]);
            % First test algorithm
            % If no FPSELECT, just use what is in CORO.FPM_ASSIGNED.

            if(nargin>1)
                if(size(FPSELECT) == CORO.FPM.size)
                    CORO.FPM_ASSIGNED = FPSELECT;
                else
                    fprintf('ERROR: FPSELECT must be the same size as the FPM');
                    return;
                end
            end
            
            BRANCH_CUT = sqrt(5); % Move to a less popular spot.

            % Initializing the Lyot field response vectors
            LPvectors = CORO.Fmatrix';
            for n=1:size(LPvectors,2)
                LPvectors(:,n) = LPvectors(:,n) * CORO.Ffp.grid_(n);
            end

            % For reference...
            NORMS = zeros(1,CORO.FPM.numel);
            for n=1:size(LPvectors,2)
                NORMS(n) = norm(LPvectors(CORO.PPMASK(:),n));
            end
            
            NORMS = NORMS(:);
            
            % Check solvability....
            V0 = LPvectors(:,CORO.FPM_ASSIGNED(:)) * exp(1i*CORO.Ffp.k*CORO.FPM.grid_(CORO.FPM_ASSIGNED(:)));
            CORO.Flyot.grid(V0).plotC(2);

            NormV0 = norm(V0);
            Usable = sum(NORMS(~CORO.FPM_ASSIGNED(:)));

            fprintf('Triangle check: \nStart=%.3g, Remainder Sum=%.3g\n',NormV0,Usable);
            if(Usable>NormV0)
                fprintf('We might be able to do this!\n');
            else
                fprintf('The triangle inequality says this is impossible.\n');
            end
            
            Nremain = sum(~CORO.FPM_ASSIGNED(:));
            ALL_NORMS = zeros(Nremain,2);
            
            while(sum(~CORO.FPM_ASSIGNED(:))>0)
                Nremain = sum(~CORO.FPM_ASSIGNED(:));
                modprint(Nremain,20);
                
                UNSELECTED = find(~CORO.FPM_ASSIGNED); % Lists numbers out of ALL pixels in the FPM.

                V0 = LPvectors(:,CORO.FPM_ASSIGNED(:)) * exp(1i*CORO.Ffp.k*CORO.FPM.grid_(CORO.FPM_ASSIGNED(:)));
                CORO.Flyot.grid(V0);
                %CORO.Flyot.plotC(4);
                %title('Partial Lyot Field (gamma 4)');

                Norm2V0 = norm(V0)^2;

                ALL_NORMS(Nremain,:) = [norm(V0),sum(NORMS(~CORO.FPM_ASSIGNED(:)))];

                if(ALL_NORMS(Nremain,1)>ALL_NORMS(Nremain,2))
                    fprintf('Triangle Inequality says IMPOSSIBLE!\n');
                    break;
                end
                
                %LPvectors = CORO.LyotContribs(~CORO.FPM_ASSIGNED); 
                PROJECT = LPvectors(CORO.PPMASK(:),UNSELECTED)'*V0(CORO.PPMASK(:));
                
                METRIC = (Norm2V0 + NORMS(UNSELECTED).^2 - 2*abs(PROJECT))./Norm2V0;
                
                %[~,best_pixel] = max(abs(PROJECT));
                %[~,best_pixel] = max(abs(PROJECT).*(abs(angle(PROJECT))>1e-3));

                [BestMetric,best_pixel] = min(METRIC);
                
                if(BestMetric>=1)
                    fprintf('Hmmmm.  All of my choices see to make things worse. Quitting.\n');
                    break;
                end
                
                NEW_PHASE = angle(-exp(-1i*BRANCH_CUT)*PROJECT(best_pixel))+BRANCH_CUT;
                fprintf('PHASE: %.3f pi. Best Metric is %f\n',NEW_PHASE/pi,BestMetric);
                
                CORO.FPM.bumpPixel1(UNSELECTED(best_pixel),NEW_PHASE/CORO.Ffp.k);
                CORO.FPM_ASSIGNED(UNSELECTED(best_pixel)) = true;

                %if(mod(sum(~CORO.FPM_ASSIGNED(:)),100)==0)
                if(true)
                    subplot(2,2,1);
                    %imagesc(CORO.FPM_ASSIGNED);sqar;axis xy;
                    %plot(ALL_NORMS','o');
                    %plot(Nremain:size(ALL_NORMS,1),ALL_NORMS(Nremain:end,:),'-');
                    semilogy(Nremain:size(ALL_NORMS,1),ALL_NORMS(Nremain:end,:),'-');
                    title('Norms');
                    subplot(2,2,2);
                    CORO.FPM.show;
                    title('FPM');
                    
                    subplot(2,2,3);
                    %CORO.Ffp_.grid(CORO.Ffp.grid())*CORO.FPM;
                    %CORO.Ffp_.plotC(4);
                    %title('post FPM field');
                    CORO.Ffp_.constant(nan).grid_(UNSELECTED)=METRIC;
                    CORO.Ffp_.show;
                    title('METRIC (and nans)');
                    
                    subplot(2,2,4);
                    CORO.Flyot * CORO;
                    %CORO.Flyot.plotC(2);
                    [x,y] = CORO.Flyot.coords;
                    imagesc(x,y,CORO.Flyot.mag2);sqar;axis xy;
                    colorbar;
                    title('Lyot Field');
                    
                    drawnow;
                end
            end
        end

        %% Taking apart the whole vector (subtraction approach).
        function CORO = algo2(CORO,FPSELECT)
            % CORO.algo1([FPSELECT]);
            % Always start from the total Lyot vector and try to reduce it
            % to zero.
            %
            % If no FPSELECT, just use what is in CORO.FPM_ASSIGNED.
            % 
            % Usage: CORO.reset(Rld>3 | Rld<0.1).algo2;

            if(nargin>1)
                if(size(FPSELECT) == CORO.FPM.size)
                    CORO.FPM_ASSIGNED = FPSELECT;
                else
                    fprintf('ERROR: FPSELECT must be the same size as the FPM');
                    return;
                end
            end
            
            % Initialize the fields for an on-axis source.
            CORO.PPtoLP(CORO.Fpp.planewave);
            %CORO.FPtoLP;
            %CORO.Flyot*CORO;
            
            BRANCH_CUT = sqrt(5); % Move to a less popular spot.

            % Initializing the Lyot field response vectors
            LPvectors = CORO.LyotContribs();

            % For reference...
            NORMS = normCols(LPvectors(CORO.PPMASK(:),:))';
            
            % Check solvability....
            % Initial Lyot Field
            V0 = LPvectors * exp(1i*CORO.Ffp.k*CORO.FPM.grid_(:));
            CORO.Flyot.grid(V0).plotC(2);

            NormV0 = norm(V0);
            Usable = sum(NORMS(~CORO.FPM_ASSIGNED(:)));

            fprintf('Triangle check: \nStart=%.3g, Remainder Sum=%.3g\n',NormV0,Usable);
            if(Usable>NormV0)
                fprintf('We might be able to do this!\n');
            else
                fprintf('The triangle inequality says this is impossible.\n');
            end
            
            Nremain = sum(~CORO.FPM_ASSIGNED(:));
            ALL_NORMS = zeros(Nremain,2);
            
            while(sum(~CORO.FPM_ASSIGNED(:))>0)
                Nremain = sum(~CORO.FPM_ASSIGNED(:));
                modprint(Nremain,20);
                
                UNSELECTED = find(~CORO.FPM_ASSIGNED); % Lists numbers out of ALL pixels in the FPM.

                V0 = LPvectors * exp(1i*CORO.Ffp.k*CORO.FPM.grid_(:));
                CORO.Flyot.grid(V0);
                %CORO.Flyot.plotC(4);
                %title('Partial Lyot Field (gamma 4)');

                Norm2V0 = norm(V0)^2;

                ALL_NORMS(Nremain,:) = [norm(V0),sum(NORMS(~CORO.FPM_ASSIGNED(:)))];

                if(ALL_NORMS(Nremain,1)>ALL_NORMS(Nremain,2))
                    fprintf('Triangle Inequality says IMPOSSIBLE!\n');
                    fprintf('RELEASING ALL PIXELS.\n');
                    
                    CORO.FPM_ASSIGNED(:) = false;
                    
                    %break;
                end
                
                %PROJECT = LPvectors(CORO.PPMASK(:),UNSELECTED)'*V0(CORO.PPMASK(:));
                PROJECT = LPvectors(CORO.PPMASK(:),:)'*V0(CORO.PPMASK(:));
                %METRIC = (Norm2V0 + NORMS(UNSELECTED).^2 - 2*abs(PROJECT))./Norm2V0;
                METRIC = (Norm2V0 + NORMS.^2 - 2*abs(PROJECT))./Norm2V0;
                [BestMetric,best_pixel] = min(METRIC(UNSELECTED));
                
                if(BestMetric>=1)
                    fprintf('Hmmmm.  All of my choices see to make things worse. Quitting.\n');
                    break;
                end
                
                %NEW_PHASE = angle(-exp(-1i*BRANCH_CUT)*PROJECT(best_pixel))+BRANCH_CUT;
                NEW_PHASE = angle(-exp(-1i*BRANCH_CUT)*PROJECT)+BRANCH_CUT;
                fprintf('PHASE: %.3f pi. Best Metric is %f\n',NEW_PHASE/pi,BestMetric);
                
                %CORO.FPM.bumpPixel1(UNSELECTED(best_pixel),NEW_PHASE/CORO.Ffp.k);
                CORO.FPM.bumpPixel1(UNSELECTED,NEW_PHASE(UNSELECTED)/CORO.Ffp.k);
                
                CORO.FPM.grid(mod(CORO.FPM.grid()-CORO.Ffp.lambda/2,CORO.Ffp.lambda))+CORO.Ffp.lambda/2;
                
                %CORO.FPM_ASSIGNED(UNSELECTED(best_pixel)) = true;
                % Allow repeated updates.
                
                %if(mod(sum(~CORO.FPM_ASSIGNED(:)),100)==0)
                if(true)
                    subplot(2,2,1);
                    %imagesc(CORO.FPM_ASSIGNED);sqar;axis xy;
                    %plot(ALL_NORMS','o');
                    %plot(Nremain:size(ALL_NORMS,1),ALL_NORMS(Nremain:end,:),'-');
                    semilogy(Nremain:size(ALL_NORMS,1),ALL_NORMS(Nremain:end,:),'-');
                    title('Norms');
                    subplot(2,2,2);
                    CORO.FPM.show;
                    title('FPM');
                    
                    subplot(2,2,3);
                    %CORO.Ffp_.grid(CORO.Ffp.grid())*CORO.FPM;
                    %CORO.Ffp_.plotC(4);
                    %title('post FPM field');
                    %CORO.Ffp_.constant(nan).grid_(UNSELECTED) = METRIC;
                    CORO.Ffp_.constant(nan).grid_(:) = METRIC(:);
                    CORO.Ffp_.show;
                    title('METRIC (and nans)');
                    
                    subplot(2,2,4);
                    CORO.Flyot * CORO;
                    %CORO.Flyot.plotC(2);
                    [x,y] = CORO.Flyot.coords;
                    imagesc(x,y,CORO.Flyot.mag2);sqar;axis xy;
                    colorbar;
                    title('Lyot Field');
                    
                    drawnow;
                end
            end
        end

        %% UNUSED
        function COEFS = antiLyotCoefs(CORO,F)
            % COEFS = antiLyotCoefs(CORO,F)
            % NOTE!!! This is currently unused.
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
        
        %%
        
        function [VMODES,s] = LyotModes(CORO,FIELD0,ThreshNum)
            % [MODES,s] = CORO.LyotModes(FIELD0,[ThreshNum]);
            % Compute the FPM modes that maximize power in the Lyot stop
            % pixels selected by CORO.PPMASK.
            % MODES are the selected svd V-modes and s is the full list of singular values.
            % ThreshNum determines how many modes to return.
            % no entry returns all of the modes.
            % ThreshNum > 1 returns that many modes.
            % ThreshNum < 1 returns those modes where s/s(1) > ThreshNum.
            
            CORO.Fpp.grid(FIELD0.grid);
            CORO.PPtoFP;
            
            MATRIX = CORO.RVmerge(CORO.Fmatrix',CORO.Ffp.grid_);
            [~,S,VMODES] = svd(MATRIX(CORO.PPMASK(:)),'econ');
            s = diag(S);
            clear S
            
            if(nargin < 3)
                return;
            end
            
            if(ThreshNum > 1) 
                VMODES = VMODES(:,1:round(ThreshNum));
            else
                VMODES = VMODES(:,s/s(1)>ThreshNum);
            end
        end
        
        function [VMODES,s] = FFPModes(CORO,FIELD0,ThreshNum)
            % [MODES,s] = CORO.FFPModes(FIELD0,[ThreshNum]);
            % Compute the FPM modes that maximize deliver power to the
            % final focal plane (FFP).
            %
            % MODES are the selected svd V-modes and s is the full list of singular values.
            % ThreshNum determines how many modes to return.
            % no entry returns all of the modes.
            % ThreshNum > 1 returns that many modes.
            % ThreshNum < 1 returns those modes where s/s(1) > ThreshNum.
            
            fprintf('Setting up the coronagraph...\n');
            CORO.Fpp.grid(FIELD0.grid);
            CORO.PPtoFP;
            
            %SELECT = (normalize(CORO.Ffp.mag2)>0.002);
            SELECT = (normalize(CORO.Ffp.mag2)>0);
            
            fprintf('Building the FP1-to-LYOT operator...\n');
            MATRIX = AOGrid.RVmerge(CORO.Fmatrix',CORO.Ffp.grid_);
            
            fprintf('Folding in the Lyot stop.\n');
            if(isempty(CORO.LYOT))
                fprintf('...using the pupil mask as the Lyot stop.\n')
                MATRIX = AOGrid.LVmerge(MATRIX,CORO.grid_(:));
            else
                fprintf('...using CORO.LYOT as the Lyot stop.\n')
                MATRIX = AOGrid.LVmerge(MATRIX,CORO.LYOT.grid_(:));
            end
            
            fprintf('Including the propagation to the final focal plane...\n');
            MATRIX = CORO.Fmatrix(SELECT,:) * MATRIX ;

            fprintf('Done building the operator.\n')

            fprintf('SVD...\n')
            %[~,S,VMODES] = svd(MATRIX,'econ');
            %[~,S,VMODES] = svd(MATRIX(SELECT,:),'econ');
            [~,S,VMODES] = svd(MATRIX,'econ');
            s = diag(S);
            clear S
            
            if(nargin < 3)
                return;
            end
            
            if(ThreshNum > 1) 
                VMODES = VMODES(:,1:round(ThreshNum));
            else
                VMODES = VMODES(:,s/s(1)>ThreshNum);
            end
        end
        
        
        %% Utilities
        function CORO = setCentroid(CORO)
            % CORO = CORO.setCentroid();
            % Set the CENTROID property for determining the optical axis.
            
            CORO.CENTROID = CORO.centroid;
        end
    end
    
    %% Static methods (utils that I want to keep in this context)
    methods(Static=true)
        
        function CTRANS = BBtrans(K,H)
            % CTRANS = AOCoronagraph.BBtrans(FREQLIST,HEIGHTS);
            % FREQLIST is the list of lambda0/lambda that you want.
            % Heights is a list of NxM heights in units of lambda0.
            % The calculation will include N pixels with M contributions each.
            
            H = sum(H,2);
            K = K(:);
            
            CTRANS = 0;
            for n=1:length(H)
                CTRANS = CTRANS + exp(1i*K*H(n));
            end
            
            CTRANS = CTRANS/length(H);
        end
        
    end    
    
end

