function CUBE = cubeReadTIFF(FILENAME)

% CUBE = cubeReadTIFF(FILENAME)

FILES = dir(FILENAME)

for n=1:length(FILES)
    try
        CUBE(:,:,n) = imread(FILES(n).name);
    catch
        break;
    end
end
