function SUM = cubeDotVec(CUBE,VEC)

% SUM = cubeDotVec(CUBE,VEC)
% 
% This is a basic datacube operation.
% SUM = SIGMA_n CUBE_n * Vec_n.  
% 
% BUGS: I don't check nothin'. 
%
% 20130719 jlcodona: Added to the NECO toolbox.

SUM = 0;

if(length(VEC) ~= size(CUBE,3))
    error('The vector has to be the same length as the 3rd dimension of the CUBE');
    return;
end

for n=1:length(VEC)
    SUM = SUM + CUBE(:,:,n)*VEC(n);
end

