function WFCUBE = cubeSubtractPiston(WFCUBE,MASK)

% WFCUBE = cubeSubtractPiston(WFCUBE,MASK)
% 

for n=1:size(WFCUBE,3)
    WAVEFRONT = WFCUBE(:,:,n);
    WAVEFRONT(MASK(:)) = WAVEFRONT(MASK(:))-mean(WAVEFRONT(MASK(:)));
    WFCUBE(:,:,n) = WAVEFRONT;
end

