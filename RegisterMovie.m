function RegisterMovie(input)
%RegisterMovie does DFT registration for the red and green channels from ThorImage. 
%Nikolas Francis 2016
tic();
for ind = 1:length(input.expname)

%Find index of experiment to be used for registration template
%regidx=find(strcmpi(input.expname,input.regexp));
%Path of images used for registration template
bb=strsplit(input.path,'/'); % KA- switced from strsep to strsplit

newpath = char(strcat(input.savepath, '/', bb{8},'/',bb{9},'/',input.expname{ind}, '/'));
%regpath1 = char(strcat(regpath, 'Experiment.xml'));

%Movies are already seperated by channel, so we only need to select 1st channel.
origpath= char(strcat(input.path, '/', input.expname{ind}, '/Experiment.xml'));
opts = get_options_from_xml(origpath);
opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};
%Load registration image sequences into memory
RedChannel = memmapfile([newpath 'redchannel.raw'], 'Format', opts.format, ...
    'Repeat', opts.numframes);
nframes = opts.numframes;
sampInd = round(linspace(1,nframes,nframes/input.winsize));
corrSeq = zeros(length(sampInd),3);
tstartCorr = tic;
%Option for finding registration template via a stable time-window, or
%first frame
if strcmpi(input.template,'stable')
    %Select red channel template from stable region of movie, ie. a set of
    %frames with high correlations
    disp('Finding "stable" segment of movie for template....');
    for k = 1:length(sampInd)-1
        %fprintf('Searching frame block %d/%d/n', k, length(sampInd)-1);
        imgSeqTmp = [];
        RedImages=[];
        for i = sampInd(k):sampInd(k+1)-1
            RedImages = cat(3,RedImages, uint8(RedChannel.Data(i).channels/256));
        end
        imgSeqTmp = reshape(RedImages, size(RedImages,1)*size(RedImages,2), ...
            size(RedImages,3));
        imgSeqTmp = double(imgSeqTmp);
        rho = corr(imgSeqTmp);
        corrSeq(k,:) = [mean(rho(:)) sampInd(k) sampInd(k+1)];
    end
    clear imgSeqTmp
    pkCorrInd = find(corrSeq(1:end-1,1) == max(corrSeq(1:end-1,1)));
    tEndCorr = toc(tstartCorr);
    disp(['Time to find stable portion of movie: ' ,num2str(tEndCorr/60), ...
        ' minutes'])
    timPtWindow = ceil((corrSeq(pkCorrInd,3) + corrSeq(pkCorrInd,2))/2);
    I=[];
    idx=floor(max([1 timPtWindow-(input.winsize/2)]):timPtWindow+(input.winsize/2));
    for i = 1:length(idx)-1
        I = cat(3,I, RedChannel.Data(idx(i)).channels);
    end
    I = squeeze(mean(I,3));
    I = padarray(I,[64 64]);
    %Spatially whiten the template, then take its Fourier transform
    template = fft2(whiten(I));
elseif strcmpi(input.template,'first')
    %Use the first image as a template
    disp('Using first frame as template....');
    I = RedChannel.Data(1).channels;
    I = padarray(I,[64 64]);
    %Spatially whiten the template, then take its Fourier transform
    template = fft2(whiten(I));
end
%Register the data
    
    %Open files for registered movies in local savepath
    bb=strsplit(input.path,'/'); % KA- changed strsep to strsplit
    newpath = [input.savepath '/' bb{8} '/' bb{9} '/' ...
        input.expname{ind} '/'];
    newpath_red = [newpath 'redchannelregistered.raw'];
    newpath_green = [newpath 'greenchannelregistered.raw'];
    outfile{1} = fopen(newpath_red, 'w');
    outfile{2} = fopen(newpath_green, 'w');
    regpath= char(strcat(input.path, '/', input.expname{ind}, '/Experiment.xml'));
    opts = get_options_from_xml(regpath);
    %Load image sequences
    opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};
    RedChannel = memmapfile([newpath 'redchannel.raw'], 'Format', opts.format, ...
        'Repeat', opts.numframes);
    GreenChannel = memmapfile([newpath 'greenchannel.raw'], 'Format', ...
        opts.format, 'Repeat', opts.numframes);
    xytrans=[];
    for frame = 1:opts.numframes
        redChanImg = RedChannel.Data(frame).channels;
        redChanImg = padarray(redChanImg,[64 64]);
        greenChanImg = GreenChannel.Data(frame).channels;
        greenChanImg = padarray(greenChanImg,[64 64]);
        %Get DFT registration coordinates. dftregistration from Mathworks
        %website.
        s = dftregistration(template,fft2(whiten(redChanImg)),input.subpixegfact);
        ty = s(3);
        tx = s(4);
        [m n]=size(redChanImg);
        %Sub-pixel registration done as a phase shift in upsampled Fourier space
        Nm = ifftshift([-fix(m/2):ceil(m/2)-1]);
        Nn = ifftshift([-fix(n/2):ceil(n/2)-1]);
        [Nn,Nm] = meshgrid(Nn,Nm);
        redChanImg = fft2(redChanImg).*exp(j*2*pi*(-ty*Nm/m-tx*Nn/n));
        redChanImg = real(ifft2(redChanImg*exp(j*s(2))));
        greenChanImg = fft2(greenChanImg).*exp(j*2*pi*(-ty*Nm/m-tx*Nn/n));
        greenChanImg = real(ifft2(greenChanImg*exp(j*s(2))));
        xytrans = [xytrans; [tx ty]];
        fwrite(outfile{1}, redChanImg(65:end-64,65:end-64), 'uint16');
        fwrite(outfile{2}, greenChanImg(65:end-64,65:end-64), 'uint16');
        %if opts.echo
            %fprintf('Registered frame %d/%d; tx = %g, ty = %g/n', frame, ...
             %   opts.numframes, tx, ty);
        %end
    end
    %Save translation coordinates
    save([newpath '/xytrans.mat'],'xytrans')
    fclose(outfile{1});
    fclose(outfile{2});
end
fprintf('Elapsed time: %g. minutes/n', toc()/60);
