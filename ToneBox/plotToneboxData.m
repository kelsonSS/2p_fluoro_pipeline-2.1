%plotToneboxData plots data for Francis et al. eNeuro. (2019)
%Nikolas Francis, 6/2019

%User name
u=getenv('username');

%code path
addpath(genpath(['C:\Users\PsiDev\Desktop\Ilan_Psignal']))

%Female cages
F = 1:2:11;

%Male cages
M = 2:2:12;

%% Habituation
phaseChoice = 1; %Training phase
trials = []; %Selected trials
Opt=0; %Optimal trials based on maximizing d'. Not available for habituation
cfact=1; %Factor for multiplying response bias
daytrials=[0 5 17]; %Select trials based on time-of-day; [A B C], A: 1 or 0 for use.
w = dir(['C:\Users\PsiDev\Desktop\Ilan_Psignal\ToneBoxData\habituation']);
Habituation=[];
parfor i = 3:length(w)
    ww=what([w(i).folder '\' w(i).name]);
    fileName = [w(i).folder '\' w(i).name '\' ww.mat{1}];
    bhvStat = toneBoxAnalysis(fileName, [],Opt,cfact,daytrials);
    Habituation=[Habituation; {bhvStat}];
        
end

%Gather data
lickOGramTarget=[];
lickOGramNonTarget=[];
LatencyH=[];
LatencyF=[];
LatencyE=[];
DoptTrialseconds=[];
DoptTrialsDays=[];
earlyCount=[];
hitCount=[];
falseAlarmCount=[];
totalTrials=[];
for i = M %1:length(Habituation)
    lickOGramTarget = [lickOGramTarget Habituation{i}.LickoGramTarget];
    LatencyH = [LatencyH Habituation{i}.lickResponseH];
    DoptTrialseconds = [DoptTrialseconds Habituation{i}.DoptTrialseconds(1,:)'];
    DoptTrialsDays = [DoptTrialsDays Habituation{i}.DoptTrialsDays(1,:)'];
    S = Habituation{i}.DoptTrialseconds(2,:);
    D = Habituation{i}.DoptTrialsDays(2,:);
    hitCount=[hitCount Habituation{i}.hitCount];
    totalTrials=[totalTrials Habituation{i}.totalTrials];
    
end
plotToneboxGroupData(phaseChoice, lickOGramTarget, lickOGramNonTarget, earlyCount, hitCount, falseAlarmCount, ...
    totalTrials, LatencyH, LatencyF, LatencyE, DoptTrialseconds, DoptTrialsDays, S, D)

%Save pdf
print('Habituation_muMales.pdf','-dpdf','-painters')
close all

%% Shaping
phaseChoice = 2;
trials = [];
Opt=0;
cfact=1;
daytrials=[0 5 17];
w = dir(['C:\Users\PsiDev\Desktop\Ilan_Psignal\ToneBoxData\shaping']);
Shaping=[];
parfor i = 3:length(w)
    ww=what([w(i).folder '\' w(i).name]);
    fileName = [w(i).folder '\' w(i).name '\' ww.mat{1}];
    bhvStat = toneBoxAnalysis(fileName, [],Opt,cfact,daytrials);
    Shaping=[Shaping; {bhvStat}];
    
end

%Gather data
lickOGramTarget=[];
lickOGramNonTarget=[];
LatencyH=[];
LatencyF=[];
LatencyE=[];
DoptTrialseconds=[];
DoptTrialsDays=[];
earlyCount=[];
hitCount=[];
falseAlarmCount=[];
totalTrials=[];
for i = 1:length(Shaping)
    lickOGramTarget = [lickOGramTarget Shaping{i}.LickoGramTarget];
    LatencyH = [LatencyH Shaping{i}.lickResponseH];
    LatencyE = [LatencyE Shaping{i}.lickResponseE];
    DoptTrialseconds = [DoptTrialseconds Shaping{i}.DoptTrialseconds(1,:)'];
    DoptTrialsDays = [DoptTrialsDays Shaping{i}.DoptTrialsDays(1,:)'];
    S = Shaping{i}.DoptTrialseconds(2,:);
    D = Shaping{i}.DoptTrialsDays(2,:);
    hitCount=[hitCount Shaping{i}.hitCount];
    earlyCount=[earlyCount Shaping{i}.earlyCount];
    totalTrials=[totalTrials Shaping{i}.totalTrials];
    
end
GroupStats=plotToneboxGroupData(2, lickOGramTarget, lickOGramNonTarget, earlyCount, hitCount, falseAlarmCount, totalTrials, ...
    LatencyH, LatencyF, LatencyE, DoptTrialseconds, DoptTrialsDays, S, D);

%Save pdf
print('Shaping_mu.pdf','-dpdf','-painters')
close all

%% Detection
phaseChoice = 3; %Training phase
trials = []; %Selected trials
Opt=0; %Optimal trials based on maximizing d'. Not available for habituation
cfact=1; %Factor for multiplying response bias
daytrials=[0 5 17]; %Select trials based on time-of-day; [A B C], A: 1 or 0 for use.
w = dir(['C:\Users\PsiDev\Desktop\Ilan_Psignal\ToneBoxData\detection']);
Detection=[];
parfor i = 3:length(w)
    ww=what([w(i).folder '\' w(i).name]);
    fileName = [w(i).folder '\' w(i).name '\' ww.mat{1}];
    disp(['Processing: ' fileName])
    bhvStat = toneBoxAnalysis(fileName, trials, Opt, cfact, daytrials);
    Detection=[Detection; {bhvStat}];
    
end

%Gather data
good = [2 4 5 6 7 8 9 10 12]
% bad = [1 3 11]

lickOGramTarget=[];
lickOGramNonTarget=[];
LatencyH=[];
LatencyF=[];
LatencyE=[];
DoptTrialsPerMinute=[];
DoptTrialsDays=[];
DoptTrialsDaysCages=[];
earlyCount=[];
hitCount=[];
falseAlarmCount=[];
totalTrials=[];
dprime=[];
for i = 1:length(Detection)
    lickOGramTarget = [lickOGramTarget Detection{i}.LickoGramTarget];
    LatencyH = [LatencyH Detection{i}.lickResponseH];
    LatencyE = [LatencyE Detection{i}.lickResponseE];
    DoptTrialsPerMinute = [DoptTrialsPerMinute Detection{i}.trialsPerMinData];
    DoptTrialsDays = [DoptTrialsDays Detection{i}.DoptTrialsDays(1,:)'];
    DoptTrialsDaysCages = [DoptTrialsDaysCages; mean(Detection{i}.DoptTrialsDaysData')];
    D = Detection{i}.DoptTrialsDays(2,:);
    hitCount=[hitCount Detection{i}.hitCount];
    earlyCount=[earlyCount Detection{i}.earlyCount];
    totalTrials=[totalTrials Detection{i}.totalTrials];
    dprime=[dprime Detection{i}.dprime];
    
end
GroupStats=plotToneboxGroupData(3, lickOGramTarget, lickOGramNonTarget, earlyCount, hitCount, falseAlarmCount, totalTrials, ...
    LatencyH, LatencyF, LatencyE, DoptTrialsDays, D, DoptTrialsDaysCages, DoptTrialsPerMinute);

%Save pdf
% print('Detection.pdf','-dpdf','-painters')
% close all

%% Discrimination
phaseChoice = 4;
trials = [];
Opt=0;
cfact=1;
daytrials=[1 17 5];
trials = [];
w = dir(['C:\Users\PsiDev\Desktop\Ilan_Psignal\ToneBoxData\discrimination']);
Discrimination=[];
for i = 4:length(w)
    ww=what([w(i).folder '\' w(i).name]);
    fileName = [w(i).folder '\' w(i).name '\' ww.mat{1}];
    bhvStat = toneBoxAnalysis(fileName, trials,Opt,cfact,daytrials);
    Discrimination=[Discrimination; {bhvStat}];
    
    %Save pdfs
    %     print(['Discrimination_' num2str(i-2) '.pdf'],'-dpdf','-painters');
    %     close all
    
end

lickOGramTarget=[];
lickOGramNonTarget=[];
LatencyH=[];
LatencyF=[];
LatencyE=[];
DoptTrialseconds=[];
DoptTrialsDays=[];
DoptTrialsecondsCages=[];
DoptTrialsDaysCages=[];
earlyCount=[];
hitCount=[];
falseAlarmCount=[];
totalTrials=[];
dprime=[];
for i = 1:length(Discrimination)
    lickOGramTarget = [lickOGramTarget Discrimination{i}.LickoGramTarget];
    LatencyH = [LatencyH Discrimination{i}.lickResponseH];
    LatencyE = [LatencyE Discrimination{i}.lickResponseE];
    LatencyF = [LatencyF Discrimination{i}.lickResponseF];
    DoptTrialseconds = [DoptTrialseconds Discrimination{i}.DoptTrialseconds(1,:)'];
    DoptTrialsDays = [DoptTrialsDays Discrimination{i}.DoptTrialsDays(1,:)'];
    DoptTrialsDaysCages = [DoptTrialsDaysCages; mean(Discrimination{i}.DoptTrialsDaysData')];
    DoptTrialsecondsCages = [DoptTrialsecondsCages; mean(Discrimination{i}.DoptTrialsecondsData')];
    S = Discrimination{i}.DoptTrialseconds(2,:);
    D = Discrimination{i}.DoptTrialsDays(2,:);
    hitCount=[hitCount Discrimination{i}.hitCount];
    earlyCount=[earlyCount Discrimination{i}.earlyCount];
    falseAlarmCount=[falseAlarmCount Discrimination{i}.falseAlarmCount];
    totalTrials=[totalTrials Discrimination{i}.totalTrials];
    dprime=[dprime Discrimination{i}.dprime];
    
end
GroupStats=plotToneboxGroupData(4, lickOGramTarget, lickOGramNonTarget, earlyCount, hitCount, falseAlarmCount, totalTrials, ...
    LatencyH, LatencyF, LatencyE, DoptTrialseconds, DoptTrialsDays, S, D, DoptTrialsDaysCages, DoptTrialsecondsCages);

%Save pdf
print('Discrimination.pdf','-dpdf','-painters')
close all