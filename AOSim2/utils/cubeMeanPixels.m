function C_ = cubeMeanPixels(CUBE)

% C_ = cubeMeanPixels(CUBE)
% Compute a mean on a per-pixel basis.
% Determine data counts by pixel ~= 0.

C_ = sum(double(CUBE),3);
N_ = sum(double(CUBE)~=0,3);

N_(N_==0) = 1; 

C_ = C_./N_;

