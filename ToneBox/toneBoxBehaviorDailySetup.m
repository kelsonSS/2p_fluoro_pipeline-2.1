%Separates tonebox detection data files into daily data files that can be
%analyzed using "toneboxDetectionDaily"
%IMPORTANT: if a detection file contains multiple stimulus types (this is
%the case for most boxes currently), each of those separate stimulus types
%are stored as separate structures within the compiled detection file ("detectionData").
%To properly execute this code, the structure corresponding to the stimulation
%type you wish to separate into daily files for analysis must be saved as a
%separate file. Ex: for 2018-2019 detection data for toneboxes 1-12, the
%compiled "detectionData" files have 3-4 structures. The structures
%containing multiple target tones and frequencies (useful for audiograms)
%have been saved as "detectFile" for each tonebox and were used to create
%the daily data files saved to each tonebox's "daily_data" folder
%Ilan Goldstein (11/2019)
box = input('Select tonebox for analysis: ');
tonebox = strcat('tonebox',num2str(box));
genPath = 'C:\Users\PsiDev\Desktop\Ilan_Psignal';                          %location of relevant psignal code and this script
dataLoc = 'D:\Devices';                                                    %parent folder containing tonebox data folders
cd(dataLoc)
%%select and load detection file to be split by day%%
[file, path] = uigetfile;                                                  %select detection file ("detectFile") to be separated for daily analysis
loadFile = fullfile(path,file);
load(loadFile)
%%list all dates associated with combined data file and find all unique dates%%
timeStamp = datetime(detectFile.timeStamp(:,2),'InputFormat','MM/dd/yyyy HH:mm:ss','Format','MM-dd-yyyy');
dateSTR = datestr(timeStamp);
for i = 1:length(dateSTR)
    dateCell{i,1} = dateSTR(i,1:11);
end
days = unique(dateCell,'stable');
dates = datetime(days,'InputFormat','dd-MMM-yyyy','Format','dd-MMM-yyyy');
%%saving daily data into separate files under "daily_data" folder%%
phaseChoice = detectFile.phaseChoice;                                      %phaseChoice, target, and nontarget should all be the same for all data
target = detectFile.target;                                                %this is why it is important to separate structures that do not have the same stimuli
nontarget = detectFile.nontarget;
n = 1;
for i = 1:length(days)
    levelVec = [];
    responseVec = [];
    toneVec = [];
    totalData = [];
    timeStamp = [];
    while strcmp(dateCell{n},days{i})
        levelVec = [levelVec; detectFile.levelVec(n)];
        responseVec = [responseVec; detectFile.responseVec(n)];
        toneVec = [toneVec; detectFile.toneVec(n)];
        totalData = [totalData; detectFile.totalData(n,:)];
        timeStamp = [timeStamp; detectFile.timeStamp(n,2)];
        n = n + 1;
        if n > length(dateCell)
            break
        end
    end
    dailyData = fullfile(dataLoc,tonebox,'daily_data',strcat(days{i},'_detection'));        %location of saved daily data files
    save(dailyData,'levelVec','toneVec','responseVec','totalData','timeStamp','phaseChoice','target','nontarget')
end
