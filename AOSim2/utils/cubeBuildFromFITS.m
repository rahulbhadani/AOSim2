function CUBE = cubeBuildFromFITS(PATTERN,PAD_SZ)

% CUBE = cubeBuildFromFITS(PATTERN,PAD_SZ)

FILES = dir(PATTERN);

if(nargin<2)
    PAD_SZ = 0;
    for n=1:length(FILES)
        HEADER = fits_read_header(FILES(n).name);
        N = double(max(HEADER.NAXIS1,HEADER.NAXIS2));
        PAD_SZ = max(PAD_SZ,N);
    end
end

CUBE = zeros([PAD_SZ,PAD_SZ,length(FILES)]);

for n=1:length(FILES)
    IMG = fits_read_image(FILES(n).name);
    IMG = padarray(IMG,PAD_SZ*[1 1]-size(IMG),'post');
    CUBE(:,:,n) = IMG;
end


