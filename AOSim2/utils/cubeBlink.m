function cubeBlink(CUBE,LIST,N,tau)

% cubeBlink(CUBE,LIST,[N],[tau])


if(nargin<2 || isempty(LIST))
    LIST = 1:size(CUBE,3);
end


if(nargin<3)
    N = 10;
end

if(nargin<4)
    tau = 0.1;
end

for nn=1:N
    for n=1:length(LIST)
        imagesc(CUBE(:,:,LIST(n)));
        sqar;
        colorbar;
        drawnow;
        
        pause(tau);
    end
end


