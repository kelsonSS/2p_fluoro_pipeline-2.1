
function [RegFname] = motionCorr3Dfxn(pthname,fname,exptVars,dftResolution,dftWinSize)

% Note: if you don't use Experiment.xml, just make sure you properly define
% dimX, dimY, and numImages
% numBytes = 2; % assumes 16-bit
segmentSize = exptVars.segmentSize; %(dimX * dimY * numBytes); %size of image in bytes
dimX = exptVars.dimX;
dimY = exptVars.dimY;
numImages = exptVars.numImages;
flybackFrames = exptVars.flybackFrames;
numZSteps = exptVars.numZSteps;
totalZplanes = exptVars.totalZplanes; %numZSteps+flybackFrames;
numVolumes = exptVars.numVolumes; %numImages/totalZplanes;

%% Load the input .raw file into a matrix
% IMG = cell(totalZplanes,1);
fullpathIMG = input.pthname fname];
for p = 1:totalZplanes
    
    IMG = zeros(dimX,dimY,numVolumes);
    
    fh = fopen(fullpathIMG); % opens the file
    fseek(fh,(p-1)*segmentSize,'bof'); % goes to the first image of this z-plane
%     fseek(fh,(p-1)*segmentSize*2,'bof'); % goes to the first image of
%     this z-plane (2-channel data)
    
    %% Obtain images of this single plan from raw file
    % loops through time, reads selected segments (planes)
    for v = 1:numVolumes
        IMG(:,:,v) = permute ( reshape((fread(fh , segmentSize/2 , 'uint16=>uint16' , 'l')) , [dimX , dimY]), [2 1] );
%         fseek(fh,(totalZplanes-1)*segmentSize*2,'cof'); % Go to desired z-plane of next volume
        fseek(fh,(totalZplanes-1)*segmentSize,'cof'); % Go to desired z-plane of next volume
%         fseek(fh,segmentSize,'cof'); % skip "C" channel
    end
    fclose(fh);
    
    %% Generate motion offsets using DFT
    offsets = calcDftOffsets(IMG,dftResolution,dftWinSize); % IMG should be MxNxFrames double
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
        tmpRegIMG = zeros(dimY+pd,dimX+pd);
        tmpRegIMG( pd/2+ty(t):pd/2+ty(t)+dimY-1 , pd/2+tx(t):pd/2+tx(t)+dimX-1) = IMG(:,:,t);
        RegIMG( : , : , t) =  tmpRegIMG(pd/2:pd/2+dimY-1 , pd/2:pd/2+dimX-1);
    end
    
    clear IMG
    
    %% Save registered .raw file
    RegFname = ['RegisteredZ' num2str(p) 'plane.raw'];
    RegFullPath = fullfile(pthname,RegFname);
    disp(['Saving Registered file in ', RegFullPath]);
    fileID = fopen(RegFullPath,'w');
    fwrite(fileID,RegIMG,'uint16','ieee-le');
    fclose(fileID);
    
    meanIMG = mean(RegIMG,3);
    save(fullfile(pthname,['RegisteredImgZ' num2str(p) '.mat']),'meanIMG','offsets','-v7.3')
%     save(fullfile(pthname,['RegisteredImgZ' num2str(p) '.mat']),'RegIMG','meanIMG','offsets','-v7.3')
    clear RegIMG
    
end

%% Can save the registered Images and offsets as .mat here if you want (probably should)
% save(fullfile(pthname,'RegisteredImg.mat'),'RegIMG','meanIMG','offsets','-v7.3')

end

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