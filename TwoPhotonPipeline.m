%%

%2 photon processing pipe-line to take raw movies from ThorImage and output fluorescence
%traces for each selected cell. It is assumed that the data will be
%processed one animal at a time, and that the stimulus delivery was done
%with Psignal. Currently, paths are based on Windows format.
% 2.0-Kelson Shilling-Scrivo 2018
% 1.0-Nikolas Francis 2016

% start initialization and file selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all; close all;

%Path to code
try 
% uncomment if you don't have github GitHub\2p_fluoro_pipeline-2.1 added and on your main path
d = dir('C:\Users\**\Documents\GitHub\2p_fluoro_pipeline-2.1');    
git_path = d(1).folder;
addpath(genpath(git_path))
cd(git_path)
catch
try
    cd('C:\Users\Kelson\Google Drive\2p_fluoro_pipeline 2.1')
    addpath(genpath('C:\Users\Kelson\Google Drive\2p_fluoro_pipeline 2.1'))
    addpath(genpath('C:\Users\Kelson\Google Drive\Psignal'))
catch
    cd('G:\My Drive\2p_fluoro_pipeline 2.1')
    addpath(genpath('G:\My Drive\2p_fluoro_pipeline 2.1'))
    addpath(genpath('G:\My Drive\Psignal'))
end

end 
%Path to Data
inpath = '\\vault3\Data\Kelson\Files to upload';
%Input variables
input=[];
%Expected frame rate form experiment. According to Thorlabs, the frame rate is
%effectivly fixed. ThorImage usually shows a frame rate + <1 Hz, which is
%why there is sometimes an extra frame at the end of a trial. However, if we
%assume that the frame estimated from the ThorImage timing.txt files to be
%the first frame of a trial is the first frame, then we can parse the movie
%into trials by taking the expected number of frames after the first frame,
%given the known trial length in time, and the expected frame rate.


input.device = 'Bscope';

input.location = 'A1';
input.age = 'p30-p60';  % !!TO-DO!!: go to Psignal and extract age and sex 
input.sex = 'm';          
input.imaging_plane = 'Layer 2_3'; % use either depth or layer ie.
                                  % 200um or Layer 2_3,
                                  % whichever is more specific
input.indicator = 'Gcamp6s';
input.excitation = 920; % in nm 


input.expectedFPS = 30;

input.emission_wavelength = 500; % GFP emission in uM


%stable = find stable window of input.winsize duration (recommended);
%first = use first frame as template
input.template = 'stable'; %'AVG_greenchannel.raw';
%# frames used to find stable window. Here default is 30 frames, ie., 1 s.
input.winsize = 1500;
%subpixel registration factor (ie. registration up to 1/input.subpixregfact)
input.subpixegfact = 10;
%neuropil (NP) subtraction factor
input.percNP = 0.7;
%expected cell dimensions for ring making
input.expectedNeuronDiamMicrons = 10;
%more dimensions, but in units of pixels
input.cellcropdim = 30/2;
input.ringthickness = 6/2;
input.NPthickness = 30/2;
input.NPgap = 3/2;
%ring smoothing factor
input.smoothfact = 30;
%plot extracted fluorescence stats
input.plotbrightness = 0;
%registration zero-pading
input.border = 1;
%path to local data storage for processed data
input.savepath = '\\Vault3\Data\Kelson\Analyzed';
%largest movie to load based on number of frames. Beyond this # the movie
%is processed in chunks for extracting fluorescence.
input.maxframechunk = 7500;


strsep = @strsplit;

% createDatalist checks for which folders don't have a corresponding
% data in ...\Analyzed and populates their paths,expnames and
% psignalfiles. matless is a logical check to ensure that all expts have
% corresponding psignalfile. NOTE: it is up to the user to determine that
% they inserted the correct psignalfile into the correct folder
[paths,expnames,psignalfiles,matless] = createDataList();
%input.animal = 'Many';
%eval([input.animal 'DataList'])
input.regexp = 'Image_0001_0001';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% start of data processing code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract red and green channels from raw images. This code will save the
%data to a subfolder within input.savepath. Subsequent registration, cell
%selection and fluourescence extraction will be done on data from within the
%input.savepath subfolder.

%%
for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Seperate red and green channels
    %ExtractRedGreenChannels(input)
     ExtractChannels(input)
end


%Register movies
for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Register combined movies
    %RegisterMovie(input)
     RegisterMovie_SingleChannel(input)
end


for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    CreateSmoothImage(input,1)
end



CreateCellDefinitionList(input) 
%%
%Click cell centers of registered movies. CellDefinitionGUI will save cell
%definitions to the directory that the images were loaded from, so be sure
%to select the images in input.savepath.
CellDefinitionGUI(input)
%%
%Generate Psignal Matrices.
for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname =          expnames(~cellfun(@isempty,expnames(:,i)),i);
    input.psignalfiles = psignalfiles(~cellfun(@isempty,psignalfiles(:,i)),i);
 
    
    
    
    for ii = 1:length(input.expname) 
        file = char(fullfile(input.path,input.expname(ii),input.psignalfiles(ii)));
    handles  = WF_getPsignalInfo(file);
    % KA- switced from strsep to strsplit
    bb = strsep(input.path, filesep);
try
    savepath = char(fullfile(input.savepath, bb{2} , input.expname(ii),'\PsignalMatrix.mat'));
    save(savepath,'handles')
catch
    savepath = char(fullfile(input.savepath, bb{end} , input.expname(ii),'\PsignalMatrix.mat'));
    disp(savepath)
    save(savepath,'handles')
end 
    end 
end



%% Extract Timing Params
for i=12:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Select current Psignal file
    input.psignalfiles = psignalfiles(~cellfun(@isempty,psignalfiles(:,i)),i);
    fprintf('analyzing %d of %d: %s \n',i,length(paths),input.path) 
    ExtractTimingParams(input,1);
end






%% Extract Fluorescence Traces 
for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Select current Psignal file
    input.psignalfiles = psignalfiles(~cellfun(@isempty,psignalfiles(:,i)),i);
    fprintf('analyzing %d of %d: %s \n',i,length(paths),input.path) 
    ExtractFluorescence(input);
end 

%% package into Experiments 
for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Select current Psignal file
    input.psignalfiles = psignalfiles(~cellfun(@isempty,psignalfiles(:,i)),i);
    fprintf('analyzing %d of %d: %s \n',i,length(paths),input.path) 
    GatherFluoroByExperiment(input);
end




%% convert files to NWB
%%  MatNWB generator initilization
 matnwb = 'C:\Users\Kelson\Documents\GitHub\matnwb';
  addpath(genpath(matnwb))
 cd(matnwb)






for i=1:length(paths)
   %Select current path
   input.path = paths{i};
   bb=strsplit(input.path,'\'); 
   input.animalID = bb{end};
   input.path = fullfile(input.savepath,input.animalID);
   input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Select current Psignal file
    input.psignalfiles = psignalfiles(~cellfun(@isempty,psignalfiles(:,i)),i);
    fprintf('Converting %d of %d: %s \n', i,length(paths),input.path)
    Mat2NWB(input)
end

  
 




%% Analysis Starts here 
 
FindFilesByExptName('FRANoise')


Data = Fluoro_to_Table_interactive('\\Vault3\Data\Kelson\Analyzed')

Data.df_by_level = getTracesAndFRA(Data) 



