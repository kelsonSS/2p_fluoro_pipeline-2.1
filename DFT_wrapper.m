

%% Edit these inputs
dftResolution = 1000;  % DFT upsampling factor

%% Load necessary experimental parameters
fullpathXML = [pthname 'Experiment.xml'];
XML = danParseXML(fullpathXML);
dimX = str2double(XML.Children(28).Attributes(16).Value);
dimY = str2double(XML.Children(28).Attributes(16).Value);
numImages = str2double(XML.Children(36).Attributes(5).Value);

%% Load the input .raw file into a matrix
fullpathIMG = [pthname fname];
fh = fopen(fullpathIMG); % or whatever the filename is
IMGtmp = permute ( reshape(fread(fh,inf,'uint16=>uint16','l'),[dimX, dimY, numImages]) , [2 1 3]);
IMG = double(IMGtmp);
fclose(fh);
clear IMGtmp

%% Generate motion offsets using DFT
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
    tmpRegIMG = zeros(dimY+pd,dimX+pd);
    tmpRegIMG( pd/2+ty(t):pd/2+ty(t)+dimY-1 , pd/2+tx(t):pd/2+tx(t)+dimX-1) = IMG(:,:,t);
    RegIMG( : , : , t) =  tmpRegIMG(pd/2:pd/2+dimY-1 , pd/2:pd/2+dimX-1);
end

%% Save registered .raw file
RegFname = 'RegisteredImage.raw';
RegFullPath = fullfile(pthname,RegFname);
disp(['Saving Registered file in ', RegFullPath]);
fileID = fopen(RegFullPath,'w');
fwrite(fileID,RegIMG,'uint16','ieee-le');
fclose(fileID);


%% Can save the registered Images and offsets as .mat here if you want (probably should)
% save('directory/filename.mat','RegIMG','offsets','-v7.3')


%% Plot side-by-side correction
% figure;
% for t = 1:size(IMG,3)
%     subplot(121)
%         imagesc(IMG(:,:,t))
%         axis('square')
%     subplot(122)
%         imagesc(RegIMG(:,:,t))
%         axis('square')
%     drawnow
% end