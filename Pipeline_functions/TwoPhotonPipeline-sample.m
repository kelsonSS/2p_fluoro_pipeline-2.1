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


% % uncomment if you don't have github GitHub\2p_fluoro_pipeline-2.1 added and on your main path
% try
% git_path = 'C:\Users\Kelson\Documents\GitHub\2p_fluoro_pipeline-2.1';
% if ~exist(git_path,'file')
%     d = dir('C:\Users\**\Documents\GitHub\2p_fluoro_pipeline-2.1');    
%     git_path = d(1).folder;
% end 
% addpath(genpath(git_path))
% cd(git_path)
% catch
% try
%     cd('C:\Users\Kelson\Google Drive\2p_fluoro_pipeline 2.1')
%     addpath(genpath('C:\Users\Kelson\Google Drive\2p_fluoro_pipeline 2.1'))
%     addpath(genpath('C:\Users\Kelson\Google Drive\Psignal'))
% catch
%     cd('G:\My Drive\2p_fluoro_pipeline 2.1')
%     addpath(genpath('G:\My Drive\2p_fluoro_pipeline 2.1'))
%     addpath(genpath('G:\My Drive\Psignal'))
% end
% 
% end 
%Input variables
input=[];

%Path to Data
input.inpath = 'Z:\\Kelson\Files to upload';

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
input.savepath = 'Z:\\Kelson\Analyzed';
%largest movie to load based on number of frames. Beyond this # the movie
%is processed in chunks for extracting fluorescence.
input.maxframechunk = 7500;


strsep = @strsplit;

% createDatalist checks for which folders don't have a corresponding
% data in ...\Analyzed and populates their paths,expnames and
% psignalfiles. matless is a logical check to ensure that all expts have
% corresponding psignalfile. NOTE: it is up to the user to determine that
% they inserted the correct psignalfile into the correct folder

input.regexp = 'Image_0001_0001.raw';

[expt_paths,psignalfiles,animalID] = createDataList(input);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for expt = 1:length(expt_paths)
    filePreparation(expt_paths{expt},input.inpath,input.savepath,psignalfiles{expt})
    
end 


% start of data processing code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract red and green channels from raw images. This code will save the
%data to a subfolder within input.savepath. Subsequent registration, cell
%selection and fluourescence extraction will be done on data from within the
%input.savepath subfolder.

%%
for expt=1:length(expt_paths)
    %Select current path
    input.expname = expt_paths{expt};
    input.animalID = animalID{expt};
    %Seperate red and green channels
    %ExtractRedGreenChannels(input)
     ExtractChannels(input)
end

    
%Register movies


for expt=1:length(expt_paths)
    %Select current path
    input.path = input.inpath;
    input.expname = expt_paths{expt};
    input.animalID = animalID{expt};
    %Register combined movies
    %RegisterMovie(input)
     RegisterMovie_SingleChannel(input)
end

% 
% for expt=1:length(expt_paths)
%     %Select current path
%     input.path = input.inpath;
%     input.expname = expt_paths{expt};
%     input.animalID = animalID(expt);
%     CreateSmoothImage(input,1)
% end

% Extract Timing Params
for expt=1:length(paths)
    %Select current path
    input.path = input.inpath;
    input.expname = expt_paths(expt);
    input.animalID = animalID(expt);
    input.psignalfiles = psignalfiles(expt);
    fprintf('analyzing %d of %d: \n',expt,length(paths)) 
    ExtractTimingParams(input,1);
end


CreateCellDefinitionList(input) 
%%
%Click cell centers of registered movies. CellDefinitionGUI will save cell
%definitions to the directory that the images were loaded from, so be sure
%to select the images in input.savepath.
CellDefinitionGUI(input)
%%
% 
badfiles = {};
for i=1:length(ExtractionPaths)
    %Select current path
    input.path = ExtractionPaths{i};
    %Select current Psignal file
    input.psignalfiles = ExtractionPsignalFiles{i};
    fprintf('analyzing %d of %d: %s \n',i,length(ExtractionPaths),input.path) 
   out =  ExtractFluorescence(input);
 %  if ~isempty(out.Errors)
  %     badfiles{end+1} = ExtractionPaths{i};
%end
end 
% end
% 




%% Extract Fluorescence Traces 

[ExtractionPaths, ExtractionPsignalFiles] = CreateCellExtractionList(input.savepath);

badfiles = {};
for expt=1:length(ExtractionPaths)
    %Select current path
    input.path = ExtractionPaths{expt};
    %Select current Psignal file
    input.psignalfiles = ExtractionPsignalFiles{expt};
    fprintf('analyzing %d of %d: %s \n',expt,length(ExtractionPaths),input.path) 
   out =  ExtractFluorescence(input,1);
   if ~isempty(out.Errors)
       badfiles{end+1} = ExtractionPaths{expt};
end
end 
%% Quality check extracted Fluorescence
for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Select current Psignal file
    input.psignalfiles = psignalfiles(~cellfun(@isempty,psignalfiles(:,i)),i);
    fprintf('QC-ing %d of %d: %s \n',i,length(paths),input.path) 
    CheckExperimentQuality(input);
end




%% package into Experiments


for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    input.expname = expnames(~cellfun(@isempty,expnames(:,i)),i);
    %Select current Psignal file
    input.psignalfiles = psignalfiles(~cellfun(@isempty,psignalfiles(:,i)),i);
    fprintf('packaging %d of %d: %s \n',i,length(paths),input.path) 
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



