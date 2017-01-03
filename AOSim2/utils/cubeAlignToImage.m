function [CUBE,SHIFTS] = cubeAlignToImage(ReferenceIMG,CUBE,VERBOSE)

% [CUBE,SHIFTS] = cubeAlignToImage(ReferenceIMG,CUBE,[VERBOSE=true])

if(nargin<3)
    VERBOSE = false;
end

SHIFTS = zeros(size(CUBE,3),2);

for n=1:size(CUBE,3)
    modprint(n,20);
    [SHIFTED,SHIFTb2a,~] = alignToImage(ReferenceIMG,CUBE(:,:,n));
    CUBE(:,:,n) = SHIFTED;
    SHIFTS(n,:) = SHIFTb2a;
end

fprintf('\n');

if(VERBOSE)
    % plot(SHIFTS(:,1),SHIFTS(:,2),'o-');sqar;drawnow;
    
    plot(sqrt(sum(SHIFTS.^2,2)),'o-');drawnow;
    
    imagesc(mean(CUBE,3));sqar;
    drawnow;
end


