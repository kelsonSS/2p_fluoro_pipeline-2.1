function RegisterMovie_SingleChannel(input)
%RegisterMovie does DFT registration for the Green and green channels from ThorImage. 

tic();
%Find index of experiment to be used for registration template
%regidx=find(strcmpi(input.expname,input.regexp));
%Path of images used for registration template
   
        
        if  isstruct(input)
            newpath = fullfile(input.savepath,input.expname) ;
              sprintf( 'registering movie %s , elapsed time: %d minutes',...
            input.expname,toc()/60)
        
        else
            newpath = input;
            input = struct('maxframechunk', 5000,'subpixegfact', 10);
        end 
        
       
disp(newpath)

% check if expt has already been registered
if ~isempty(dir(fullfile(newpath, 'greenchannelregistered*')))
    return
end 

% if ~ isempty(dir(fullfile(newpath, 'greenchannelregistered.raw')))
%     newfile = dir(fullfile(newpath, 'greenchannelregistered.raw'));
%     oldfile = dir(fullfile(newpath, 'greenchannel.raw'));
%     if newfile.bytes == oldfile.bytes
%         return
%     end
% end

 if isempty(dir(fullfile(newpath, 'greenchannel.raw')))
    return
 end 

%Movies are already seperated by channel, so we only need to select 1st channel.
origpath= fullfile(newpath,'Experiment.xml');
opts = get_options_from_xml(origpath);
opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};
%Load registration image sequences into memory


%% Load the input .raw file into a matrix
FullPathIMG = fullfile(newpath, 'greenchannel.raw');

fh = fopen(FullPathIMG); 
numFrames_total =  FindRawImgSize(fh,[opts.dimX opts.dimY]);
ItemsPerImage = opts.dimX* opts.dimY ;
chunkSize =  input.maxframechunk; % frames 
chunks = ceil(numFrames_total/chunkSize); 
fclose(fh);
tstartCorr = tic;


%if opts.totalZplanes > 1
%    Register3dMovie(FullPathIMG,opts,newpath)
%    return
%end 

% 2d Data registration
for chunk_count = 1:chunks
    %% register data 
     
    % open file and move pointer to next chunk
    fh = fopen(FullPathIMG); 
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

txty = calcDftOffsets(IMG,input.subpixegfact);
txty= round(txty);

% winsize = 99; 
% nframes = size(IMG,3);
% sampInd = round(linspace(winsize+1,nframes-winsize,10));
% corrSeq = zeros(length(sampInd),3);
% 
% 
% 
% if ~exist('imTemplate','var')
%     disp('Finding "stable" segment of movie....');
%     tstartCorr = tic;
%     for k = 1:length(sampInd)
%         lwin = max(1, sampInd(k)-winsize);
%         rwin = min(sampInd(k)+winsize, nframes);
%         imgSeqTmp = reshape(IMG(:,:,lwin:rwin), size(IMG,1)*size(IMG,2),size(IMG(:,:,lwin:rwin),3));
%         rho = corr(imgSeqTmp);
%         corrSeq(k,:) = [mean(rho(:)) lwin rwin];
%     end %K
%     clear imgSeqTmp
%     
%     pkCorrInd = find(corrSeq(1:end-1,1) == max(corrSeq(1:end-1,1)));
%     tEndCorr = toc(tstartCorr);
%     fprintf(' Time to find stable portion of movie %d seconds \n',round(tEndCorr))
%     
%     
%    
%     timPtWindow = (corrSeq(pkCorrInd,3) + corrSeq(pkCorrInd,2))/2;
%     lwin = max(1, timPtWindow-winsize);
%     rwin = min(timPtWindow+winsize, nframes);
%     
%     I = (mean(IMG(:,:,lwin:rwin),3));
%   
%         
%     fixed = (I - min(I(:)))./range(I(:));
%     imTemplate = fft2(fixed);
%     save(fullfile(newpath,'AvgImg.m'),'fixed')
% end
% if isempty(gcp('nocreate'))
%     parpool('local');
% end
%  %% perform DFT registration & get motion correction offsets
% disp('')
% disp('Finding motion correction coordinates....');
% tStartMotionOffsets = tic;
% parfor j = 1 :  nframes
%     % using Fourier transformation of images for registration
%     error  = dftregistration(imTemplate,fft2(IMG(:,:,j)),10);
%     ty(j) = error(3);
%     tx(j) = error(4);
% end
% 
% txty = [ty' tx'];
% txty = round(txty);
% 
% telapsedMotionOffset = toc(tStartMotionOffsets);
% disp('')
% disp(['     Time Elapsed for Motion Offsets was: ', num2str(telapsedMotionOffset), ' seconds'])
%  

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
    RegFullPath = fullfile(newpath, RegFname); 
    %RegFullPath = fullfile('C:\Users\KanoldLab\Desktop\Kelson\registration_temp', RegFname);
    fprintf('Saving Registered file %d of %d in %s \n',chunk_count,chunks, RegFullPath);
    if chunk_count == 1
        fileID = fopen(RegFullPath,'w+');
    else % append if in batch mode  
        fileID = fopen(RegFullPath,'a+');
    end 
    fwrite(fileID,RegIMG,'uint16','ieee-le');
    fclose(fileID);
    fprintf('saved! loading next file');
   

%% Plot side-by-side correction
%  figure;
%  for ii =1:1000
% imagesc(imfuse(fixed,RegIMG(:,:,ii)))
% pause(.01)
% end
    
    
    
clear IMG
clear RegIMG
clear tmpRegIMG
clear txty
clear tx
clear ty

%  if chunk_count == chunks
%          RegFname = 'greenchannelregistered.raw';
%          TransferFullPath = [newpath RegFname];
%         movefile(RegFullPath,TransferFullPath)
%         deletefile(RegFullPath)
%  end


end
clear imTemplate



%% Can save the registered Images and offsets as .mat here if you want (probably should)
% save('directory/filename.mat','RegIMG','tx','ty','-v7.3')



    
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

%fclose(fh);
 

fprintf('Elapsed time: %g. minutes \n', toc()/60);
end 


function Register3dMovie(FullPathIMG,opts,newpath)


dimX = opts.dimX;
dimY = opts.dimY;
segmentSize = dimX * dimY * 2;    %size of image in bytes

flybackFrames = opts.flybackFrames;
numZSteps = opts.numZSteps;
totalZplanes = opts.totalZplanes; %numZSteps+flybackFrames;
fh = fopen(FullPathIMG);
numImages = FindRawImgSize(fh);
fclose(fh);
numVolumes = floor(numImages / totalZplanes);

for p = 1:totalZplanes
    
    IMG = zeros(dimX,dimY,numVolumes);
    
    fh = fopen(FullPathIMG); % opens the file
    fseek(fh,(p-1)*segmentSize,'bof'); % goes to the first image of this z-plane
%     fseek(fh,(p-1)*segmentSize*2,'bof'); % goes to the first image of
%     this z-plane (2-channel data)
    
    %% Obtain images of this single plane from raw file
    % loops through time, reads selected segments (planes)
    for v = 1:numVolumes
        IMG(:,:,v) = permute ( reshape((fread(fh , segmentSize/2 , 'uint16=>uint16' , 'l')) , [dimX , dimY]), [2 1] );
%         fseek(fh,(totalZplanes-1)*segmentSize*2,'cof'); % Go to desired z-plane of next volume
        fseek(fh,(totalZplanes-1)*segmentSize,'cof'); % Go to desired z-plane of next volume
%         fseek(fh,segmentSize,'cof'); % skip "C" channel
    end
    fclose(fh);
    
    %% Generate motion offsets using DFT
    dftResolution = 10;
   offsets = calcDftOffsets(IMG,dftResolution); % IMG should be MxNxFrames double
    ty = round(offsets(:,1));
    tx = round(offsets(:,2));
    
    %% Apply offsets
    dimY = size(IMG,1);
    dimX = size(IMG,2);
    nframes = size(IMG,3);
    pdTmp = round(max(dimY,dimX)/3); % 1/3 of IMG size for padding
    pd = pdTmp + mod(pdTmp,2); % Make padding amount an even number
    RegIMG = zeros( dimY , dimX , nframes);
    for t = 1:nframes
        
   try
        tmpRegIMG = zeros(dimY+pd,dimX+pd);
        tmpRegIMG( pd/2+ty(t):pd/2+ty(t)+dimY-1 , pd/2+tx(t):pd/2+tx(t)+dimX-1) = IMG(:,:,t);
        RegIMG( : , : , t) =  tmpRegIMG(pd/2:pd/2+dimY-1 , pd/2:pd/2+dimX-1);
   catch
        RegIMG(:,:,t) = IMG(:,:,t);
        warning('Problem registering frame %d',t)
    end
         
   end
    
    clear IMG
    
    %% Save registered .raw file
    RegFname = ['GreenChannelRegisteredZ' num2str(p) 'plane.raw'];
    RegFullPath = fullfile(newpath,RegFname);
    disp(['Saving Registered file in ', RegFullPath]);
    fileID = fopen(RegFullPath,'w');
    fwrite(fileID,RegIMG,'uint16','ieee-le');
    try 
        fclose(fileID);
    catch
        fclose all
    end
    
    meanIMG = mean(RegIMG,3);
    save(fullfile(newpath,['RegisteredImgZ' num2str(p) '.mat']),'meanIMG','offsets','-v7.3')
%     save(fullfile(newpath,['RegisteredImgZ' num2str(p) '.mat']),'RegIMG','meanIMG','offsets','-v7.3')
    clear RegIMG
    
end

%% Can save the registered Images and offsets as .mat here if you want (probably should)
% save(fullfile(newpath,'RegisteredImg.mat'),'RegIMG','meanIMG','offsets','-v7.3')


%% Plot side-by-side correction -- edited for saving a GIF. Have to keep IMG and RegIMG in workspace
% gifFilename = 'motionCorr3Dexample.gif';
% figure;
% for t = 1:size(IMG,3)
%     subplot(121)
%         imagesc(IMG(:,:,t))
%         colormap('gray')
%         axis('square')
%     subplot(122)
%         imagesc(RegIMG(:,:,t))
%         colormap('gray')
%         axis('square')
%     drawnow
%     
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     
%     if t == 1
%         imwrite(imind,cm,gifFilename,'gif', 'Loopcount',inf,'DelayTime',0.1);
%     else
%         imwrite(imind,cm,gifFilename,'gif','WriteMode','append','DelayTime',0.1);
%     end
% end

end 