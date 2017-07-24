function F = propagateFarField(F,Z)
  
% AOField.propagateFarField(Z)
% Propagate the AOField assuming the Fresnel Scale sqrt(lambda*Z) is bigger
% than the support of the field.  This uses Fourier optics to go from a
% planar field to a modulated spherical wave.  The AOField sampling is
% automatically changed to a larger value.  
% 
% Z: propagation distance in m. (May be positive or negative.)

POWER0 = F.totalPower;

F.spacing(F.lambda./F.extent*Z);

% [X,Y] = F.COORDS;
% SAG = F.dsphere(Z,X,Y);

F.grid(F.fft);

F.scaleTotalPower(POWER0);

% F.grid(F.fft.*(exp(-1i*F.k*SAG)/sqrt(prod(F.size))));
% F.grid(F.fft./(sqrt(prod(F.size))));

