function C_ = cubeMean(CUBE,N)

% C_ = cubeMean(CUBE,N)
% compute a downsampled cube using means of sequences.

NC = size(CUBE,3);

if(nargin<2)
    NORM = sum(CUBE~=0,3);
    NORM(NORM==0) = 1;
    C_ = sum(CUBE,3)./NORM;
    return;
else
    N_ = floor(NC/N);
    C_ = zeros(size(CUBE(:,:,1:N_)));
    
    for n=1:N_
        SELECT = N*(n-1) + (1:N)
        C = CUBE(:,:,SELECT);
        
        NORM = sum(C~=0,3);
        NORM(NORM==0) = 1;
        
        C_(:,:,n) = sum(C,3)./NORM;
    end
end
