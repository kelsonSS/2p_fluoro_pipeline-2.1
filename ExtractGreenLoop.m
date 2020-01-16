
%2 photon processing pipe-line to take raw movies from ThorImage and output fluorescence
%traces for each selected cell. It is assumed that the data will be
%processed one animal at a time, and that the stimulus delivery was done
%with Psignal. Currently, paths are based on Windows format.
%Nikolas Francis 2016

clear

%Path to code
addpath(genpath('/Users/test/Desktop/Kevin'))

%Input variables
input=[];
%Expected frame rate form experiment. According to Thorlabs, the frame rate is
%effectivly fixed. ThorImage usually shows a frame rate + <1 Hz, which is
%why there is sometimes an extra frame at the end of a trial. However, if we
%assume that the frame estimated from the ThorImage timing.txt files to be
%the first frame of a trial is the first frame, then we can parse the movie
%into trials by taking the expected number of frames after the first frame,
%given the known trial length in time, and the expected frame rate.
input.expectedFPS = 30;
%stable = find stable window of input.winsize duration (recommended); 
%first = use first frame as template
input.template = 'stable';
%# frames used to find stable window. Here default is 30 frames, ie., 1 s.
input.winsize = 1000;
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
input.plotbrightness = 1;
%registration zero-pading
input.border = 20;
%path to local data storage for processed data
input.savepath = '/Volumes/Vault2Data/Kevin/Tone_Rearing/BoxControl';
%largest movie to load based on number of frames. Beyond this # the movie
%is processed in chunks for extracting fluorescence.
input.maxframechunk = 15000;

strsep = @strsplit;
%Load pre-constructed datalist for the animal
input.animal = 'Many';
eval([input.animal 'DataList'])
input.regexp = 'Image_0001_0001';

%Extract red and green channels from raw images. This code will save the
%data to a subfolder within input.savepath. Subsequent registration, cell
%selection and fluourescence extraction will be done on data from within the
%input.savepath subfolder.
for i=1:length(paths)
    %Select current path
    input.path = paths{i};
    %Seperate red and green channels
    ExtractRedChannel(input)
end