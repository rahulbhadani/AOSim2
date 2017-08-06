function F = propagateFarField(F,deltaZ,newExtent,newSpacing)

% AOField.propagateFarField(Z,[newExtent],[newSpacing])
% Propagate the AOField assuming the Fresnel Scale sqrt(lambda*Z) is bigger
% than the support of the field.  This uses Fourier optics to go from a
% planar field to a modulated spherical wave.  The AOField sampling is
% automatically changed to a larger value.
%
% Z: propagation distance in m. (May be positive or negative.)
% The new grid params default to the current values.
% The power should be conserved before the newExtent trimming.

POWER0 = F.totalPower;

% This preserves the initial grid params.
Fin = F.copy;
Fin.interpolate_method = 'cubic';

Fin.FFTSize = Fin.size;
Fin.spacing(F.lambda./F.extent*deltaZ);

if(F.direction>0)
    Fin.grid(Fin.clearCache.fft);
else
    Fin.grid(Fin.clearCache.ifft);
end

Fin.scaleTotalPower(POWER0);

% F output grid is different.
if(nargin>2)
    F.extent(newExtent);
end
if(nargin>3)
    F.spacing(newSpacing);
end

F.zero + Fin;  % This cubic interpolates to the Fout grid.  No spherical wave yet.

[X,Y] = F.COORDS;
SAG = F.dsphere(deltaZ,X,Y);
F.grid(F.grid.*(-1i*F.direction*exp(-1i*F.direction*F.k*SAG)));

F.z = F.z + F.direction*deltaZ;

