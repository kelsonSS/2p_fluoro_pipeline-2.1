%% Main script to analyze volumetric movies (dF/F specifically)
% Zac Bowen
%% Note: to batch analyze multiple datasets you could easily make a cell
%   full of directory names and then put the entire code in a for-loop.
%  Example:

% datadir =  {'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031020_367n_100um20st_FRA\' % 1
%             'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031020_367r_100um20st_FRA\' % 2
%             'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031020_ia09_100um20st_FRAinj\' % 3
%             'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031120_352ll_100um20st_FRA\' % 4
%             'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031120_352ll_100um20st_FRA_diffxy\'}; % 5
%         
% for fileNum = 1:length(datadir)
%     clearvars -except datadir fileNum
%     exptDir = datadir{fileNum};
%     cd(exptDir)
%   
%     --CODE--
% 
% end

%%
close all
clear

%% Define path with input .raw file & Experiment.xml file
exptDir = 'C:\Users\Zac\Documents\raw volume data\031020\367n\100um20st_FRA\';
cd(exptDir)

%% For truncating files to save RAM if testing
numImagesOverwrite = [];  %[] for all images
if ~isempty(numImagesOverwrite); fprintf('Using %.0f frames\n',numImagesOverwrite); end

%% Load necessary experimental parameters
% Note: if you don't use Experiment.xml, just make sure you properly define
% dimX, dimY, and numImages
fullpathXML = [exptDir 'Experiment.xml'];
XML = danParseXML(fullpathXML);
exptVars = xmlVersionCheck(XML,numImagesOverwrite);

%% Use Kelson's trick to find numImages if the xml file is wrong
fname = 'Image_0001_0001.raw';  % Name of input .raw file
fh = fopen(fname);
exptVars.numImages = FindRawImgSize(fh);
exptVars.numVolumes = floor(exptVars.numImages/exptVars.totalZplanes);
fclose(fh);

%% Loop through z-planes and compute DFF
imageRate = 30/exptVars.totalZplanes; % In Hertz
for z = 1:exptVars.totalZplanes
    fprintf('\n\n\n------Starting Z Plane #%.0f------\n',z)
    %% Calculating dff from Registered images
    inputImgFile = ['RegisteredZ' num2str(z) 'plane.raw'];  % Name of input .raw file
    exptName = ['Z' num2str(z)];
    
    % computeDFF below will look for a specific CellDef file naming scheme (['CellDefFile_' exptName '.mat'])
    winSizeSeconds = 10; percentBaselineSub = 50; neuropilSubPercent = 70;
    analyzedOutputFile = computeDFF(exptDir,exptName,inputImgFile,exptVars,winSizeSeconds,percentBaselineSub,neuropilSubPercent);
     
end
