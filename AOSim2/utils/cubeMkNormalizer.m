function NNN = cubeMkNormalizer(CUBE)

% NNN = cubeMkNormalizer(CUBE)
% 

AVG = mean(CUBE,3);
NNN = median(AVG(:))./AVG;

