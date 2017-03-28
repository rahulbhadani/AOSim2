classdef AOShackHartmann < AOWFS
    % AOShackHartmann < AOWFS
    %
    % This class implements a Shack-Hartmann WFS.
    % The grid is the SPOTS image and is the size of the last AOField.
    % The coordinates of the lenslets are defined when calling .mkLenslets.
    % I plan on adding methods that make it easy to set up a Fried geometry
    % configuration.
    
    properties(GetAccess='public',SetAccess='public')
        FL_lenslets = 0;
        qscale = 1;  % This is for quiver plots.
    end
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        LENSLETS = []; % This is an AOScreen lens.
        DiamSubAp = 0;  % This is for the ad hoc version.
        
        MASKED = []; % True if the subap is masked by the pupil.
        maskThresh = 0.5;
        
        Xsubap = [];  % Coords of the lenslet subaps.
        Ysubap = [];
        
        XSLOPES = []; % Slopes of the wavefront at the subaps.  In arcsecs. (horiz)
        YSLOPES = []; % (vertical)
        
        XSLOPES0 = [];  % If defined, these are subtracted.
        YSLOPES0 = [];
    end
    
    methods
        % Constructor
        function SHWFS = AOShackHartmann(arg1)
            SHWFS = SHWFS@AOWFS(arg1);  % Just define the grid.
        end
        
        function WFS = sense(WFS,F)
            % WFS = sense(WFS,F)
            % This uses a lenslet array if it is defined, otherwise it uses
            % a fake method that looks at field patches.
            
            if(isempty(WFS.LENSLETS))
                WFS.senseFake1(F);
            else
                error('The real LENSLET version of the AOShackHartmann still has to be finished.')
            end
        end
        
        function WFS = defineSubAps(WFS,X,Y,diam)
            % WFS = defineSubAps(X,Y,[diam])
            
            if(nargin<3)
                WFS.DiamSubAp = diam;
            end
            
            WFS.Xsubap = X; % Note that I leave this in whatever shape is given me.
            WFS.Ysubap = Y;
        end
        
        function WFS = defineSubApsByAperture(WFS,A,subap_spacing,COORDS0)
            % SHWFS.defineSubApsByAperture(A,subap_spacing,[COORDS0=[0 0]])
            % I make sure one of the subaps is on COORDS0.
            
            if(nargin<4)
                COORDS0 = fliplr(A.Offset);
            end
            
            BB = A.BBox;
            xmin = BB(1,2);
            xmax = BB(2,2);
            ymin = BB(1,1);
            ymax = BB(2,1);
            
            % This is overkill, but I trim the locations back.
            EXTENT = A.extent;
            Nx =  ceil(EXTENT(1)/subap_spacing);
            Ny =  ceil(EXTENT(2)/subap_spacing);
            
            WFS.DiamSubAp = subap_spacing;
            
            % Overkill and trim.
            x = (-Nx:Nx)*subap_spacing + COORDS0(1);
            y = (-Ny:Ny)*subap_spacing + COORDS0(2);
            
            x(x<xmin | x>xmax) = [];
            y(y<ymin | y>ymax) = [];
            
            [X,Y] = meshgrid(x,y);
            
            WFS.Xsubap = X(:);
            WFS.Ysubap = Y(:);
            
            WFS.setMasked(A,0.1);
        end
        
        
        function WFS = setBias(WFS)
            WFS.XSLOPES0 = WFS.XSLOPES;
            WFS.YSLOPES0 = WFS.YSLOPES;
        end
        
        function WFS = setMasked(WFS,A,threshold)
            % WFS = setMasked(WFS,A,[threshold=0.5])
            % Set the supap masks based on an AOAperture mask.
            
            if(nargin>2)
                WFS.maskThresh = threshold;
            end

            Acopy = AOSegment(A);  % Note that we don't want this to be an AOAperture.
            Acopy.padBy(ceil(WFS.DiamSubAp/A.dx));
            Acopy.grid(Acopy.convolve(fspecial('disk',WFS.DiamSubAp/A.dx/2)));
            
            X = WFS.Xsubap;
            Y = WFS.Ysubap;
            
            aper = gather(Acopy.interpGrid(X,Y));
            
            WFS.MASKED = (aper < WFS.maskThresh);
        end
        
        function n = nSubAps(WFS)
            n = sum(~WFS.MASKED(:));
        end
        
        
        function [Sx,Sy] = slopeArrays(WFS)
            if(isempty(WFS.XSLOPES0) || isempty(WFS.YSLOPES0))
                Sx = WFS.XSLOPES;
                Sy = WFS.YSLOPES;
            else
                Sx = WFS.XSLOPES - WFS.XSLOPES0;
                Sy = WFS.YSLOPES - WFS.YSLOPES0;
            end
        end
        
        function WFS = plotSubAps(WFS,linespec)
            % WFS.plotSubAps([linespec])
            
            if(nargin<2)
                linespec = 'b-';
            end
            
            N = numel(WFS.Xsubap);
            
            for n=1:N
                hold on;
                drawCircles(WFS.DiamSubAp/2,[WFS.Xsubap(n),WFS.Ysubap(n)],1,linespec);
                if(WFS.MASKED(n))
                    hold on;
                    plot(WFS.Xsubap(n),WFS.Ysubap(n),'rx');
                end
            end
            hold off;
        end
        
        function [X,Y] = subAps(WFS)
            % [X,Y] = WFS.subAps()
            % Return the coords of the unmasked subaps.
            
            X = []; Y = [];
            if(isempty(WFS.Xsubap) || isempty(WFS.Ysubap))
                return;
            end
            
            if(isempty(WFS.MASKED))
                X = WFS.Xsubap; Y = WFS.Ysubap;
            else
                X = WFS.Xsubap(~WFS.MASKED);
                Y = WFS.Ysubap(~WFS.MASKED);
            end
        end
        
        function WFS = senseFake1(WFS,F)
            % WFS = senseFake1(WFS,F)
            % This is a fake SH.  The real version should use a lenslet
            % array if it is defined.
            
            [FCUBE,dx] = WFS.extractFieldPatches(F);
            FCUBE(isnan(FCUBE)) = 0;
            
            dwf1 = angle(squeeze(mean(mean(FCUBE(2:end,:,:).*conj(FCUBE(1:end-1,:,:))))))*F.lambda/2/pi;
            dwf2 = angle(squeeze(mean(mean(FCUBE(:,2:end,:).*conj(FCUBE(:,1:end-1,:))))))*F.lambda/2/pi;
        
            WFS.XSLOPES = 206265 * dwf2 / dx;
            WFS.YSLOPES = 206265 * dwf1 / dx;
            
            if(length(WFS.XSLOPES) ~= length(WFS.XSLOPES0))
                WFS.XSLOPES0 = zeros(size(WFS.XSLOPES));
                WFS.YSLOPES0 = zeros(size(WFS.YSLOPES));
            end
        end        
        
        function [FCUBE,dx] = extractFieldPatches(WFS,F)
            % [FCUBE,dx] = WFS.extractFieldPatches(F);
            % This is for the FAKE ShackHartmann mode.
            % This only returns patches for the unmasked regions.
            % dx is the spacing of the field patches.

            d = WFS.DiamSubAp;
            Nf = ceil(d/F.dx/2)*2+1;
            x = linspace(-d/2,d/2,Nf);
            dx = x(2)-x(1);
            [Xf,Yf] = meshgrid(x);

            [Xwfs,Ywfs] = WFS.subAps;
            
            FCUBE = zeros([size(Xf),numel(Xwfs)]);

            for n=1:numel(Xwfs)
                FCUBE(:,:,n) = gather(F.interpGrid(Xf+Xwfs(n),Yf+Ywfs(n)));
            end
        end        
        
        function SLOPES = slopes(WFS)
        % SLOPES = WFS.slopes();
        % Returns a vector of slopes from the last sense measurement.
        
            SLOPES = [WFS.XSLOPES(:)-WFS.XSLOPES0(:); WFS.YSLOPES(:)-WFS.YSLOPES0(:)];
        end
                
        function WFS = quiver(WFS,overplot,linespec)
            % WFS.quiver(overplot,linespec)
            
            if(nargin<2)
                overplot = true;
            end
            
            if(nargin<3)
                linespec = 'b-';
            end
            
            [X,Y] = WFS.subAps;
            [Sx,Sy] = WFS.slopeArrays;
            
            if(overplot)
                hold on;
            end
            quiver(X(:),Y(:),Sx(:),Sy(:),linespec);
        end

    end
end

