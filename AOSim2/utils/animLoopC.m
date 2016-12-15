function animLoopC(dcube,nloops,GAMMA,RNG,secsDelay)

% animLoopC(dcube,nloops,GAMMA,RNG,secsDelay)
%
% 20050718: JLC: added support for complex data cubes.
% 20160817 JLC: Made better control for complex.

if(nargin<2)
    nloops = 1
end

if(nargin<3)
    GAMMA = 2;
end

if(nargin<4)
    MIN = min(abs(dcube(:)));
    MAX = max(abs(dcube(:)));
else
    MIN = RNG(1);
    MAX = RNG(2);
end

if(MIN <= MAX)
    MIN = MIN - 1;
    MAX = MAX + 1;
end

% set(gcf,'DoubleBuffer','on');
% set(gcf,'BackingStore','on');

%colormap(gray(256).^.5);

[NX NY NT] = size(dcube);

for iloop=1:nloops
    for iframe=1:NT
        plotComplex(squeeze(dcube(:,:,iframe)),GAMMA);
        %daspect([1 1 1]);
        axis xy
        title(sprintf('frame %d',iframe));
        drawnow;
%         if(nargin>4)
%             pause(secsDelay);
%         end
    end
end
