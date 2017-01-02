classdef AOCoronagraph < AOSegment
    % AOCoronagraph class.
    % Works basically like an AOSegment.
    %
    % 20151224: JLCodona.  UA.SO.CAAO.AOSim2
    
    properties
        REMAPPED_APERTURE    = [];
        APODIZER    = [];
        FPM         = [];
        LYOT        = [];
       
        CENTROID = [];
        
        verbose     = false; % print debugging info.
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
        
    end
end

