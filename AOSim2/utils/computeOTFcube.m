function OCUBE = computeOTFcube(CUBE,CEN)

% OCUBE = computeOTFcube(CUBE,CEN)    

OCUBE = zeros(size(CUBE));

for n=1:size(CUBE,3)
    if(nargin<2)
        OCUBE(:,:,n) = computeOTF(demean(CUBE(:,:,n)));
    else
        OCUBE(:,:,n) = computeOTF(demean(CUBE(:,:,n)),CEN);
    end
end

