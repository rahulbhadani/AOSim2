function C_ = cubeMeanPixelsNaN(CUBE)

% C_ = cubeMeanPixelsNaN(CUBE)
% Compute a mean on a per-pixel basis.
% Determine data counts by pixel ~= nan.

C_ = sum(double(CUBE),3);
N_ = sum(~isnan(CUBE),3);

N_(N_==0) = 1; 
C_ = C_./N_;

