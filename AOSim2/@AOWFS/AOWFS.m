classdef AOWFS < AOGrid
    % AOWFS Class
    % This is a simple non-physical (abstract) WFS model.  
    % Use it as a base class for making specific WFS models.
    % 20090421 JLCodona.  Part of AOSim2.
    % 20170322 JLC.  Rewriting this based on my original philosophy.
    %
    % Note that Wavefront Sensors are very different and their object
    % models should also be.  In general, the AOWFS base class doesn't
    % force any model on the sensor, just abstract plotting, a general
    % AOWFS.sense(AOField); method, and some magic static methods.  
    % 
    % This reboot will break older AO simulations, but it will make for a
    % better and more useful future.  I am not the one who took this code
    % off-course, but I was resonsible for allowing it to persist.  
    % I have learned a lot about saying "no" tp code contributitors over the years. 
    % I could rant on for some time about how understanding object models
    % is not that hard, and how breaking the model is a sign that you are
    % not qualified for bigger things.  But I fear the world is full of
    % people breaking things they don't understand.
    %
    % Virtual nod to Angry Linus.  Hey, I'm right there with you buddy!

    properties(GetAccess='public',SetAccess='public')
    end
    
    properties(GetAccess = 'public', SetAccess = 'protected')
	end
    
    methods
        % Constructor
        function WFS = AOWFS(arg1)
            % AOWFS Constructor.
            % The meaning of the grid is undefined.  
            % It should be used for something that makes sense for the WFS,
            % but it is specific to the sensor.  
            % The only uniformly required methods are the constructor and
            % WFS.sense(FIELD).  The rest are specific to the WFS type.
            
            obj = WFS@AOGrid(arg1); % This does nothing by default.
        end
        
        function WFS = sense(WFS,F)
            % WFS.sense(FIELD);
            % This returns itself, so an implementation should provide
            % methods that return the measurements and the derived results.

            error('error: The WFS.sense(FIELD) method needs to be implemented.');
        end
        
        function WFS = initBias(WFS,A)
            % WFS.initBias(A)
            % This method is WFS-specific and should take a reference
            % measurement to set the output bias.  
            % Because it includes an assumption of use into a general
            % method, I have decided to deprecate it.  
            
            error('error: The WFS.initBias method is deprecated.');
        end
        
        function WFS = setBias(WFS)
            % WFS.setBias()
            % This method sets the last WFS-specific measurement as the zero bias.
            % Subsequent measurements will include the bias.
            % Implementors should think about how this should most
            % logically work for the specific WFS type.

            error('error: The WFS.setBias method needs to be implemented.');
        end
    end
    
    
    %% static methods
    methods(Static=true)
        function TipTilt = senseTipTilt(F)
            % TipTilt = senseTipTilt(F)
            % This is a very basic PSF centroid WFS.
            % The result is given in arcsecs.
            
            TipTilt = gather(centroid(abs(F.fft).^2,F.FAXIS_PIXEL).*F.dk/F.k*206265);
            
        end
        
    end % static methods
    
end

