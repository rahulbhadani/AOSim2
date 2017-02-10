function F = propagate2(F,Z,MaxTheta)
  
% AOField.propagate2(dz,MaxTheta)
% 
% dz: propagation distance in m.
% MaxTheta: max included angles in arcsecs.
% 
% 20150424 JLCodona

  NFRESNEL = 10;
  LATERAL = 1/4;

  if(nargin<3)
      MaxTheta = sqrt(F.lambda/Z) * NFRESNEL;
      MaxTheta = min(MaxTheta,min(F.extent)*LATERAL/Z)*206265;
  end
  
  if(nargin < 2)
    error('requires at least 2 arguments: propagate(AOField,distance).');
  end
  
  % pad the field by some number of Fresnel scales.
  Rf = sqrt(abs(Z)*F.lambda);
  %fprintf('DEBUG: The Fresnel scale of this jump is %g m.\n',Rf);
  
  SZ = F.size;
  
  n1 = 1:SZ(1);
  n1 = fftshift(n1);
  n1 = n1-n1(1);

  n2 = 1:SZ(2);
  n2 = fftshift(n2);
  n2 = n2-n2(1);
  
  DK = 2*pi./F.extent;
  
  [K1,K2] = meshgrid(DK(1)*n1,DK(2)*n2);
  Knyquist = min(DK .* F.size / 2);
  
  Kcutoff = min(0.85*Knyquist,F.k*MaxTheta);
  
  % Note that this is corner-centered by construction.
  PROPAGATOR = exp(-1i*Z*F.dsphere(F.k,K1,K2)); 
  
  if(isempty(F.cache.LPF))
      KR2 = K1.^2 + K2.^2;
      %F.cache.LPF = (abs(K1)<Kcutoff).*(abs(K2)<Kcutoff);
      %F.cache.LPF = exp(-((K1.^2+K2.^2)/Kcutoff^2).^(4/2)); % Super Gaussian.
      F.cache.LPF = exp(-((KR2)/Kcutoff^2).^(4/2)); % Super Gaussian.
  end

  PROPAGATOR = PROPAGATOR .* F.cache.LPF;
  
  %F.grid(ifftshift(ifft2(PROPAGATOR.*fft2(fftshift(F.grid))))); 
  F.grid((ifft2(PROPAGATOR.*fft2((F.grid))))); 
  F.z = F.z - Z;
