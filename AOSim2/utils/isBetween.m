function yesno = isBetween(X,a,b)
% yesno = isBetween(X,a,b)
% This is an interior interval check to see if X \in (a,b).
% a and b can be listed in any order.
% X can be an array, resulting in a boolean array output.

LIMS = sort([a b]);
yesno = (X>LIMS(1)) & (X<LIMS(2));

