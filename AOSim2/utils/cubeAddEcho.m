function CUBE = cubeAddEcho(CUBE,OFFSET,FACTOR)

% CUBE = cubeAddEcho(CUBE,OFFSET,FACTOR)
% 
% Do a circshift on each frame and add an echo.
% 
% BUGS: I use circshift. 
%
% 20130609 jlcodona: Added to the NECO toolbox.

for n=1:size(CUBE,3)
    CUBE(:,:,n) = CUBE(:,:,n) ...
        + FACTOR*circshift(CUBE(:,:,n),OFFSET);
end
