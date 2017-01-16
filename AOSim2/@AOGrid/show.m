function AOG = show(AOG,RANGE)
% AOGrid.show([RANGE]);  
% Make an image of the AOGrid.

if(nargin<2)
    AOG.plotC;
    title(AOG.describe,'FontSize',14);
    daspect([1 1 1]);
else
    AOG.plotC(RANGE);
    title(AOG.describe,'FontSize',14);
    daspect([1 1 1]);
end

% AOG.center;
% 
% if(AOG.isX)
%     AOG.plotC;
%     title(AOG.describe,'FontSize',14);
%     daspect([1 1 1]);
% else
%     [x,y] = AOG.coords;
%     imagesc(x,y,AOG.ndex,[-3 0]);
%     axis square;
%     axis xy;
%     title(AOG.describe,'FontSize',14);
%     colorbar;
%     daspect([1 1 1]);
% end

end

