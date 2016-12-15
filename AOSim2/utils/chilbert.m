function HF = chilbert(F)
% HF = chilbert(F)
% Hilbert transform for complex inputs.
% 20161206 JLCodona

HF = hilbert(real(F));
HF = HF + 1i*hilbert(imag(F));


