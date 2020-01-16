function CreateSmoothImage(input,debug) 

% This function creates an smoothed average image of the movie with 
% clearly visible cells for use in collecting cell definitions.
% Specifically this function smooths the image, removing PMT shot noise from 
% the movie. Then it takes a max image of the resulting movie and saves it.
% Kelson Shilling-Scrivo 2019 
tic()
if ~exist('debug','var')
    debug = 0 ;
end 
for expnum = 1:length(input.expname)
    sprintf( 'registering movie\n %s \n elapsed time: %d minutes \n',...
    fullfile(input.path,input.expname{expnum},floor(toc()/60) ))
    bb=strsplit(input.path,'\'); % KA- switced from strsep to strsplit
    
    newpath = char(strcat(input.savepath, '\', bb{end} , '\' , input.expname{expnum}, '\' )) ;
    % check if expt has already been registered
 if ~debug
    if ~ isempty(dir([newpath 'AvgImageSmooth.tif']))
             continue
     end
 end 
     opts = get_options_from_xml(fullfile(newpath,'Experiment.xml'));
    
    fh = fopen(fullfile(newpath,'greenchannelregistered.raw'));
    IMG = fread(fh,opts.dimX* opts.dimY *1000,'uint16=>uint16');
    IMG = reshape(IMG,opts.dimX, opts.dimY,[]);
    fclose(fh);
    
    %  bin values together to remove PMT noise 
    bin_size = 10;
    r = mod(size(IMG,3), bin_size); 
    IMG = IMG(:,:,1:end-r);
    
    IMG = reshape(IMG,opts.dimX, opts.dimY,bin_size,[]);
    smooth_movie = squeeze(sum(IMG,3)) /bin_size;
      % max filter and hm filter
    smooth_img = squeeze(max(smooth_movie,[],3));
    smooth_img = hmfilter(smooth_img);
    clear IMG
    clear smooth_movie
    % ploting
    h = figure;imagesc(smooth_img);colormap gray;axis square;axis off
    saveas(h,fullfile(newpath,'AvgImageSmooth.tif'))
    save(fullfile(newpath,'smooth_img.mat'),'smooth_img')
    clear smooth_img;
   
    
end 

function [avgfilteredImage filteredImages] = hmfilter(I)
filteredImages=zeros(size(I));
for i = 1:size(I,3)
    I_temp = I(:,:,i);
    I_temp = im2double(I_temp);
    I_temp = log(1 + I_temp);
    M = 2*size(I_temp,1);
    N = 2*size(I_temp,2);
    sigma = 1.5; 
    [X, Y] = meshgrid(1:N,1:M);
    centerX = ceil(N/2);
    centerY = ceil(M/2);
    gaussianNumerator = (X - centerX).^2 + (Y - centerY).^2;
    H = exp(-gaussianNumerator./(2*sigma.^2));
    H=1-H;
    H = fftshift(H);
    Ir = padarray(I_temp,[ceil(size(I_temp,1)/2) ceil(size(I_temp,2)/2)],'symmetric');
    If = fft2(Ir, M, N);
    Iout = real(ifft2(H.*If));
    Iout = Iout(ceil(size(I_temp,1)/2)+1:size(Iout,1)-ceil(size(I_temp,1)/2), ...
        ceil(size(I_temp,2)/2)+1:(size(Iout,2)-ceil(size(I_temp,2)/2)));
    filteredImages(:,:,i) = exp(Iout) - 1;
end
avgfilteredImage = squeeze(mean(filteredImages,3));



