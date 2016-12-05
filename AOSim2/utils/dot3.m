function SUM = dot3(CUBE,vec)

% SUM = dot3(CUBE,vec)

SUM = 0;
for n=1:length(vec)
    SUM = SUM + CUBE(:,:,n)*vec(n);
end
