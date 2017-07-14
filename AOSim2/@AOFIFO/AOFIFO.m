%% AOFIFO.  Support class for AOSim2.
% 20090425 JLCodona.

classdef AOFIFO < handle
    %AOFIFO A general delay line class for AOSim2.
    %  20090425 JLCodona.
    
    properties
        pipeline = {};
        
        name = '<FIFO>';
        contains = [];    % Set this to a class name. It will be tested in push.
        % Empty means don't test.
    end
    
    methods
        function FIFO = AOFIFO(Length,classtype)
            FIFO.pipeline = cell(1,Length);
            if(nargin>1)
                FIFO.contains = classtype;
            else
                FIFO.contains = [];
            end
        end
        
        function OUTPUT = push(FIFO,INPUT)
            % This is the main operation of the FIFO.
            
            if(~isempty(FIFO.contains))
                if(~isa(INPUT,FIFO.contains))
                    fprintf('AOFIFO: WARNING: INPUT is a %s not a %s.\n',...
                        class(INPUT),FIFO.contains);
                end
            end
            
            OUTPUT = FIFO.pipeline{end};
            FIFO.pipeline = circshift(FIFO.pipeline,[0 1]);
            FIFO.pipeline{1} = INPUT;
        end
        
        function FIFO = flush(FIFO)
            FIFO.pipeline = cell(size(FIFO.pipeline));
        end
    end
end

