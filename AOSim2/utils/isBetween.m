function yesno = isBetween(X,a,b,inclusive)
% yesno = isBetween(X,a,b,[inclusive=false])
% This is an interior interval check to see if X \in (a,b).
% a and b can be listed in any order.
% X can be an array, resulting in a boolean array output.
% The optional inclusize flag determines if endpoints count.

if(nargin<4)
    inclusive = false;
else
    inclusive = logical(inclusive);
end

LIMS = sort([a b]);

if(inclusive)
    yesno = (X>=LIMS(1)) & (X<=LIMS(2));
else
    yesno = (X>LIMS(1)) & (X<LIMS(2));
end
