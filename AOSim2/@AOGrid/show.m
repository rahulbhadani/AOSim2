function AOG = show(AOG,RANGE)
% AOGrid.show([RANGE]);
% Make an image of the AOGrid.

if(nargin<2)
    AOG.plotC;
    title(AOG.describe,'FontSize',10);
    daspect([1 1 1]);
else
    AOG.plotC(RANGE);
    title(AOG.describe,'FontSize',10);
    daspect([1 1 1]);
end


end

