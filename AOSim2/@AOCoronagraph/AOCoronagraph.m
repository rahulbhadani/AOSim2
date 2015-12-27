classdef AOCoronagraph < AOSegment
    % AOCoronagraph class.
    % Works basically like an AOSegment.
    %
    % 20151224: JLCodona.  UA.SO.CAAO.AOSim2
    
    properties
        APERTURE    = [];
        APODIZER    = [];
        FPM         = [];
        LYOT        = [];
        
        verbose     = false; % print debugging info.
    end
    
    methods
        
        % Constructor
        function CORO = AOCoronagraph(varargin)
            CORO = CORO@AOSegment(varargin);
        end
        
        
    end
end

