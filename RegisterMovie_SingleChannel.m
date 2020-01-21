function RegisterMovie_SingleChannel(input)
%RegisterMovie does DFT registration for the Green and green channels from ThorImage. 

tic();
%Find index of experiment to be used for registration template
%regidx=find(strcmpi(input.expname,input.regexp));
%Path of images used for registration template
for expnum = 1:length(input.expname)
    sprintf( 'registering movie %s , elapsed time: %d minutes',...
            input.expname{expnum},toc()/60)
bb=strsplit(input.path,'\'); % KA- switced from strsep to strsplit

newpath = char(strcat(input.savepath, '\', bb{end} , '\' , input.expname{expnum}, '\' )) ;
disp(newpath)

% check if expt has already been registered
if ~ isempty(dir([newpath 'greenchannelregistered.raw']))
    newfile = dir([newpath 'greenchannelregistered.raw']);
    oldfile = dir([newpath 'greenchannel.raw']);
    if newfile.bytes == oldfile.bytes
        continue
    end
end
   
%Movies are already seperated by channel, so we only need to select 1st channel.
origpath= char(strcat(input.path, '\', input.expname(expnum), '\Experiment.xml'));
opts = get_options_from_xml(origpath);
opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};
%Load registration image sequences into memory


%% Load the input .raw file into a matrix
fullpathIMG = [newpath 'greenchannel.raw'];

fh = fopen(fullpathIMG); 
numFrames_total =  FindRawImgSize(fh,[opts.dimX opts.dimY]);
ItemsPerImage = opts.dimX* opts.dimY ;
chunkSize =  input.maxframechunk; % frames 
chunks = ceil(numFrames_total/chunkSize); 
fclose(fh);

tstartCorr = tic;
for chunk_count = 1:chunks
    %% register data 
     
    % open file and move pointer to next chunk
    fh = fopen(fullpathIMG); 
    start_idx = (chunk_count-1) * ItemsPerImage * chunkSize * 2 %bytes  ;
    fseek(fh,start_idx,'bof')
    % check to ensure we don't go over IMG size on last img
    if chunk_count == chunks 
        IMG = fread(fh,inf,'uint16');
    else 
        IMG =fread(fh,ItemsPerImage * chunkSize,'uint16');
    end 
    fclose(fh);


    IMG =  permute( reshape(IMG,opts.dimX, opts.dimY, []),[2 1 3]); 
    

%% Generate motion offsets using DFT
%% FIND SEGMENT OF MOVIE WITH HIGH CORR VALUES
winsize = 99; 
nframes = size(IMG,3);
sampInd = round(linspace(winsize+1,nframes-winsize,10));
corrSeq = zeros(length(sampInd),3);



if ~exist('imTemplate','var')
    disp('Finding "stable" segment of movie....');
    tstartCorr = tic;
    for k = 1:length(sampInd)
        lwin = max(1, sampInd(k)-winsize);
        rwin = min(sampInd(k)+winsize, nframes);
        imgSeqTmp = reshape(IMG(:,:,lwin:rwin), size(IMG,1)*size(IMG,2),size(IMG(:,:,lwin:rwin),3));
        rho = corr(imgSeqTmp);
        corrSeq(k,:) = [mean(rho(:)) lwin rwin];
    end %K
    clear imgSeqTmp
    
    pkCorrInd = find(corrSeq(1:end-1,1) == max(corrSeq(1:end-1,1)));
    tEndCorr = toc(tstartCorr);
    fprintf(' Time to find stable portion of movie %d seconds \n',round(tEndCorr))
    
    
   
    timPtWindow = (corrSeq(pkCorrInd,3) + corrSeq(pkCorrInd,2))/2;
    I = (mean(IMG(:,:,timPtWindow-winsize:timPtWindow+winsize),3));
    fixed = (I - min(I(:)))./range(I(:));
    imTemplate = fft2(whiten(fixed));
    save(fullfile(newpath,'AvgImg.m'),'fixed')
end
if isempty(gcp('nocreate'))
    parpool('local');
end
 %% perform DFT registration & get motion correction offsets
disp('')
disp('Finding motion correction coordinates....');
tStartMotionOffsets = tic;
parfor j = 1 :  nframes
    % using Fourier transformation of images for registration
    error  = dftregistration(imTemplate,fft2(whiten(IMG(:,:,j))),10);
    ty(j) = error(3);
    tx(j) = error(4);
end

txty = [ty' tx'];
txty = round(txty);

telapsedMotionOffset = toc(tStartMotionOffsets);
disp('')
disp(['     Time Elapsed for Motion Offsets was: ', num2str(telapsedMotionOffset), ' seconds'])
 

%% Apply offsets
dimY = size(IMG,1);
dimX = size(IMG,2);
nframes = size(IMG,3);
pdTmp = round(max(dimY,dimX)/3); % 1/3 of IMG size for padding
pd = pdTmp + mod(pdTmp,2); % Make padding amount an even number

RegIMG = zeros( dimY , dimX , nframes );

%% register data 
% check to ensure we don't go over IMG size on last img
tRegistration= tic; 


for t = 1:nframes
    tmpRegIMG = zeros(dimY+pd,dimX+pd);
    try
        tmpRegIMG( pd/2+txty(t,1):pd/2+txty(t,1)+dimY-1 ,...
            pd/2+txty(t,2):pd/2+txty(t,2)+dimX-1) = IMG(:,:,t);
        RegIMG( : , : , t) =  tmpRegIMG(pd/2:pd/2+dimY-1 , pd/2:pd/2+dimX-1);
    catch % in case there is a large jump in the frame larger than padding
        %tmpRegIMG(pd/2:pd/2+dimY-1,pd/2:pd/2+dimX-1) = IMG(:,:,t)
        RegIMG(:,:,t) = IMG(:,:,t);
        warning('Problem registering frame %d',t)
    end
end

telapsedRegistration = toc(tRegistration);
fprintf('Registration took %d seconds \n',ceil(telapsedRegistration))
    
    %% Save registered .raw file
    RegFname = 'greenchannelregistered.raw';
    RegFullPath = fullfile('C:\Users\KanoldLab\Desktop\Kelson\registration_temp', RegFname);
    fprintf('Saving Registered file %d of %d in %s \n',chunk_count,chunks, RegFullPath);
    if chunk_count == 1
        fileID = fopen(RegFullPath,'w+');
    else % append if in batch mode  
        fileID = fopen(RegFullPath,'a+');
    end 
    fwrite(fileID,RegIMG,'uint16','ieee-le');
    fclose(fileID);
    fprintf('saved! loading next file');
   
clear IMG
clear RegIMG
clear tmpRegIMG
clear txty
clear tx
clear ty

 if chunk_count == chunks
        RegFname = 'greenchannelregistered.raw';
        TransferFullPath = [newpath RegFname];
        movefile(RegFullPath,TransferFullPath)
        deletefile(RegFullPath)
 end


end




%% Can save the registered Images and offsets as .mat here if you want (probably should)
% save('directory/filename.mat','RegIMG','tx','ty','-v7.3')


%% Plot side-by-side correction
%  figure;
%  for ii =1:1000
% imagesc(imfuse(fixed,RegIMG(:,:,ii)))
% pause(.01)
% end
    
%%

% 
% 
% % check extraction was successful
% if size(IMG,3) ~= opts.numframes
%     
%     warnpath = fopen(fullfile(input.savepath,bb{end},input.expname{expnum},'Warnings.txt'),...
%                      'w+');
%     fprintf(warnpath,' Warning: File length does not match expected length \n');
%     
%     fclose(warnpath);
% end 
%% cleanup

fclose(fh);

end 

fprintf('Elapsed time: %g. minutes \n', toc()/60);
end 
