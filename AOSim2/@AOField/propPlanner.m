function [SPACING,EXTENT,MAX_ANGLE] = propPlanner(F,D,ROI,RANGE,DO_IT)
% [SPACING,EXTENT,MAX_ANGLE] = propPlanner(F,D,ROI,RANGE,[DO_IT])
% 
% Make suggestions for AOField size to propagate a distance RANGE.
% D is the extent of the illuminated initial field.
% ROI is the size of the final field over which we want decent accuracy.
% RANGE is the planned propagation distance.
% Note: MAX_ANGLE is in arcsecs.
% SPACING is the largest value.  Make it smaller.
% EXTENT can be larger, but it may be wasteful.
% Use MAX_ANGLE in the propagation methods.
% The grid size will be at least NxN where N is EXTENT./SPACING.

MAX_RAD = abs((D+ROI)/RANGE);
MAX_ANGLE = MAX_RAD * 206265;

SPACING = F.lambda/MAX_RAD;
EXTENT = max(abs(D),abs(ROI))+abs(MAX_RAD*RANGE);

if(nargin>4 && DO_IT)
    F.spacing(SPACING);
    F.resize(ceil(EXTENT/SPACING));
    F.clearCache;
end

