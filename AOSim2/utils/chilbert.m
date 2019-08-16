function H = chilbert(F)
% H = chilbert(F)
% Hilbert transform for complex-valued functions.
% 20161207: JLCodona, University of Arizona

H = hilbert(real(F)) + 1i * hilbert(imag(F));


