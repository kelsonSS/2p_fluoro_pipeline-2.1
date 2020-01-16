function [bhvStat bhvOutput] = toneBoxAnalysis(fileName, trials, Opt, cfact, useDayTrials)
%% toneBoxAnalysis analyses and plots data from individual ToneBoxes
% Nikolas Francis, 6/2019

% fileName: optional file path
% trials: analysis trials
% cfact: multiplication factor for optimal trial criterion
% useDayTrials: select hours for analyzed trials

if nargin < 5
    useDayTrials=0;
end

if nargin < 4
    cfact = 1;
end

if nargin < 3
    Opt = 0;
end

figure('position',[14 649 1884 311])
activitywindow=100;

%% Initialize variables and load data
totalData = []; % lick traces: time x trial
responseVec = []; % H: hit, M: miss, F: false alarm, R: correct rejection, E: early
toneVec = []; % presented tone indexes or frequencies
levelVec = []; % presented tone levels
timeStamp = []; % time stamp for each trial
H = []; % binary vector indicating hits
E = []; % binary vector indicating early licks
F= []; % binary vector indicating false alarms
Frqs = num2str(round(1.0*2.^[0:.5:5.5]',1)); % frequencies used in the toneBox
Lvls = num2str([-30:5:0]'); % levels used in the toneBox

% Load data
if isempty(fileName)
    [fileName filePath] = uigetfile('*.mat','Load tonebox data');
    fullFile = [filePath fileName];
    load(fullFile);
else
    load(fileName);
end
% make sure frequency and level presentations are in matrix not cell
if iscell(toneVec)
    toneVec = cell2mat(toneVec);
end
if iscell(levelVec)
    levelVec = cell2mat(levelVec);
end

% No tones in habituation and longer trial
if phaseChoice == 1
    Tidx = [0 10];
    toneVec=num2str(zeros(size(responseVec)));
    toneVec= toneVec(find(~isspace(toneVec)));
else
    Tidx = [0 4];
end

% Selected trials
TotalOrigTrials = length(responseVec);
if nargin<2 || isempty(trials)
    T = [1 length(responseVec)];
else
    T = trials;
end
if length(T) == 1
    T = [1 T];
end
T=sort(T);
if max(T) > length(responseVec)
    T(end) = length(responseVec);
end
if T(1) < 1
    T(1) = 1;
end

% Assemble early, hit, miss, false alarm and correct reject trials vectors
E=zeros(length(responseVec),1);
H=zeros(length(responseVec),1);
F=zeros(length(responseVec),1);
for i = 1:length(responseVec)
    if strcmpi(responseVec(i),'E')
        E(i)=1;
    elseif strcmpi(responseVec(i),'H')
        H(i)=1;
    elseif strcmpi(responseVec(i),'F')
        F(i)=1;
    end
end

%% Corrections for initial method of labeling data (No longer needed for data collected after 12/08/18
if isempty(toneVec) || isempty(totalData)
    close(gcf)
    bhvStat = [];
    bhvOutput = [];
    return
end
[toneVec toneVecN levelVec levelVecN] = correctToneLevelFreq(Frqs, Lvls, toneVec, levelVec);
if isempty(toneVecN)
    close(gcf)
    bhvStat = [];
    bhvOutput = [];
    return
end
%Reconstruct totalData until data collection is fixed in phaseChoice==4, i.e., tone frequency discrimination
if phaseChoice ==4
%     totalData=[];
    iT=0;
    iNT=0;
    %change target and nontarget frequencies to Frqs index identifier used in toneVec
    for i = 1:length(Frqs)
        if strcmpi(regexprep(Frqs(i,:),' ',''),num2str(target))
            target = i;
        elseif strcmpi(regexprep(Frqs(i,:),' ',''),num2str(nontarget))
            nontarget = i;
        end
    end
    for i = 1:length(toneVec)
        if strcmpi(toneVec(i),num2str(target))
            iT = iT+1;
%             totalData(i,:) =  totalDataTarget(iT,:);
            totalDataTarget(iT,:) = totalData(i,:);
        elseif  strcmpi( toneVec(i),num2str(nontarget))
            iNT = iNT+1;
%             totalData(i,:) =  totalDataNonTarget(iNT,:);
            totalDataNonTarget(iNT,:) = totalData(i,:);
        end
    end
end

%% Plot trial by trial licking
subplot(1,8,2:4)
cla
hold on
area(1:length(H),100*movmean(H,activitywindow),'facecolor','b','facealpha',1,'edgecolor','none');
xlim([1 length(responseVec)+1])
aa=axis;
ylim([aa(3) aa(4)+5])
set(gca,'fontsize',8)
if phaseChoice == 4
    area(1:length(E),100*movmean(E,activitywindow),'facecolor','g','facealpha',1,'edgecolor','none');
    area(1:length(F),100*movmean(F,activitywindow),'facecolor','r','facealpha',1,'edgecolor','none');
    h=legend('T','NT','E','autoupdate','off');
elseif phaseChoice ~= 1
    area(1:length(E),100*movmean(E,activitywindow),'facecolor','g','facealpha',1,'edgecolor','none');
    h=legend('T','E','autoupdate','off');
else phaseChoice == 1
    h=legend('T','autoupdate','off');
end
legend boxoff
title('Trialwise Responses')
xlabel('Trials')
set(gca,'box','off')

%% Select trials
responseVec = responseVec(T(1):T(2));
totalData = totalData(T(1):T(2),:);
toneVec = toneVec(T(1):T(2));
levelVec = levelVec(T(1):T(2));
timeStamp = timeStamp(T(1):T(2),:);
H = H(T(1):T(2));
E = E(T(1):T(2));
F= F(T(1):T(2));

% Day/Night trials
if useDayTrials(1)
    formatIn = 'MM/dd/yyyy HH:mm:ss';
    Times=datetime;
    for i = 1:length(timeStamp)
        DateString = timeStamp{i,2};
        Times(i) = datetime(DateString,'InputFormat',formatIn);
    end
    Daytrials=[];
    for i = 1:length(Times)
        if useDayTrials(3) > useDayTrials(2)
            if Times(i).Hour > useDayTrials(2) && Times(i).Hour <= useDayTrials(3)
                Daytrials = [Daytrials; i];
            end
        else
            if Times(i).Hour <= useDayTrials(3) || Times(i).Hour > useDayTrials(2)
                Daytrials = [Daytrials; i];
            end
        end
    end
    responseVec = responseVec(Daytrials);
    totalData = totalData(Daytrials,:);
    toneVec = toneVec(Daytrials);
    levelVec = levelVec(Daytrials);
    toneVecN = toneVecN(Daytrials);
    levelVecN = levelVecN(Daytrials);
    timeStamp = timeStamp(Daytrials,:);
    H = H(Daytrials,:);
    E = E(Daytrials,:);
    F= F(Daytrials,:);
end

% Optimal trials
if Opt && phaseChoice ~=1
    Data.toneVec=toneVec;
    Data.responseVec=responseVec;
    Data.totalData=totalData;
    Data.target = target;
    Data.H=H;
    Data.E=E;
    Data.F=F;
    if phaseChoice==4
        Data.nontarget = nontarget;
    end
    [optTrials c] = OptimalToneboxTrials(Data, activitywindow, cfact);
    
    % Plot criterion
    subplot(1,8,2:4)
    plot([1 TotalOrigTrials],[c c].*100,'color','k','linewidth',2)
    if phaseChoice == 4
        h=legend('T','NT','E','c','autoupdate','off');
    else
        h=legend('T','E','c','autoupdate','off');
    end
    legend boxoff
    if T(2)-T(1) ~= TotalOrigTrials-1
        area(T,repmat(aa(4)*.5,[1 length(T)]),'edgecolor','none','facecolor','k','facealpha',.2)
        text(mean(T),(aa(4)*.5)+.75,'Analysis trials','fontsize',7,'HorizontalAlignment','center')
    end
else
    optTrials = 1:length(responseVec);
end
responseVec = responseVec(optTrials);
totalData = totalData(optTrials,:);
timeStamp = timeStamp(optTrials,:);
H = H(optTrials,:);
E = E(optTrials,:);
F= F(optTrials,:);
hitCount=sum(H);
earlyCount=sum(E);
falseAlarmCount=sum(F);
totalTrials = length(responseVec);
try
    toneVec = toneVec(optTrials);
    levelVec = levelVec(optTrials);
    toneVecN = toneVecN(optTrials);
    levelVecN = levelVecN(optTrials);
catch
    toneVec=zeros(1,length(optTrials));
    levelVec=zeros(1,length(optTrials));
    toneVecN =zeros(1,length(optTrials));
    levelVecN = zeros(1,length(optTrials));
    
end

%% Plot the lick-o-gram, ie., the average lick response
subplot(1,8,1)
cla
if phaseChoice ~=4
    L = 100*mean(totalData);
    LL = smooth(L,10);
    TrialTime = linspace(0,Tidx(2),length(LL));
    bhvStat.LickoGramTarget = LL;
    plot(TrialTime,LL,'b','linewidth',2);
else
    L = 100*mean(totalDataTarget);
    LL = smooth(L,10);
    TrialTime = linspace(0,Tidx(2),length(LL));
    bhvStat.LickoGramTarget = LL;
    plot(TrialTime,LL,'b','linewidth',2);
    hold on
    L = 100*mean(totalDataNonTarget);
    LL = smooth(L,10);
    bhvStat.LickoGramNonTarget = LL;
    plot(TrialTime,LL,'r','linewidth',2);
end
xlim(Tidx);
hold on
title('Lick-o-gram')
set(gca,'fontsize',8)
aa=axis;
ylim([aa(3) aa(4)])
if phaseChoice ~= 4
    h=legend('T','AutoUpdate','off','location','northwest');
    legend boxoff
else
    h=legend('T','NT','AutoUpdate','off','location','northwest');
    legend boxoff
end
if sum(toneVec) & phaseChoice ~=1
    %Box around target
    y1 = [0 aa(4)]*.5;
    x1 = [1 1];
    y2 = [aa(4) aa(4)]*.5;
    x2 = [2 2];
    area([x1 x2],[y1 y2],'facecolor','k','facealpha',.25,'edgecolor','none');
end
ylabel('Percent (%)')
xlabel('Time(s)')
set(gca,'box','off')

%% Plot response latency, aka first lick
subplot(1,8,5)
cla

nbins = size(totalData,2);
lickResponseE = zeros(1,nbins);
lickResponseH = zeros(1,nbins);
lickResponseF = zeros(1,nbins);
firstLick=[];
for i = 1:size(totalData,1)
    firstLick=find(totalData(i,:)>0,1);
    if H(i)
        lickResponseH(1,firstLick) = lickResponseH(1,firstLick) + 1;
    elseif F(i)
        lickResponseF(1,firstLick) = lickResponseF(1,firstLick) + 1;
    elseif E(i)
        lickResponseE(1,firstLick) = lickResponseE(1,firstLick) + 1;
    end
    firstLick=[];
end
LickSum = (lickResponseE+lickResponseH+lickResponseF);
t=linspace(0,Tidx(2),length(LickSum));

if phaseChoice ~=1
    L = 100*lickResponseE/sum(LickSum);
    LL = smooth(L,10);
    Eidx = find(t<=1);
    bhvStat.lickResponseE = L(Eidx)';
    area(t(Eidx),LL(Eidx),'facecolor','g','facealpha',.5,'edgecolor','none');
    hold on
    L = 100*lickResponseH/sum(LickSum);
    LL = smooth(L,10);
    idx = find(t>=1);
    idx = [min(idx)-1 idx];
    bhvStat.lickResponseH = L(idx)';
    area(t(idx),LL(idx),'facecolor','b','facealpha',.5,'edgecolor','none');
    L = 100*lickResponseF/sum(LickSum);
    LL = smooth(L,10);
    bhvStat.lickResponseF = L(idx)';
    area(t(idx),LL(idx),'facecolor','r','facealpha',.5,'edgecolor','none');
else
    L = 100*lickResponseH/sum(LickSum);
    LL = smooth(L,10);
    idx = find(t>=0);
    bhvStat.lickResponseH = L(idx)';
    area(t(idx),LL(idx),'facecolor','b','facealpha',.5,'edgecolor','none');
    hold on
end
aa=axis;
ylim([aa(3) aa(4)*1.1])
xlabel('Time(s)')
xlim([0 Tidx(end)]);
hold on
title('Response Latency')
set(gca,'fontsize',8)
aa=axis;
if sum(toneVec) & phaseChoice ~=1
    %Box around target
    y1 = [0 aa(4)]*.5;
    x1 = t([idx(1) idx(1)]);
    y2 = [aa(4) aa(4)]*.5;
    x2 = [2 2];
    area([x1 x2],[y1 y2],'facecolor','k','facealpha',.25,'edgecolor','none');
end
set(gca,'box','off')

%% Plot tone frequency response spectrum (i.e., behavioral frequency response function/area)
nospec=0;
if phaseChoice  ~= 1 & length(unique(toneVec))>1
    subplot(1,8,6)
    cla
    freqs = 1:length(Frqs);
    levels = [0 1:length(Lvls)-1];
    R=[];
    for i = 1:length(freqs)
        if sum(levelVecN)>0
            for ii= 1:length(levels)
                idx = find(toneVecN==freqs(i) & levelVecN==levels(ii));
                R(ii,i)=sum(responseVec(idx) =='H' | responseVec(idx)=='F')./length(toneVec);
            end
        else
            idx = find(toneVecN==freqs(i));
            R(1,i)=sum(responseVec(idx) =='H' | responseVec(idx)=='F')./length(toneVec);
        end
    end
    
    bhvStat.freqs = freqs;
    bhvStat.levels = levels;
    bhvStat.Frqs = Frqs;
    bhvStat.Lvls = Lvls;
    bhvStat.freqResponse = 100*R;
    
    if ~sum(levelVecN)
        plot(1:length(freqs),100*R,'k');
        hold on
        for i = 1:length(freqs)
            if freqs(i)==target
                plot(i,100*R(i),'bs','markerface','b')
            elseif exist('nontarget') && freqs(i)==nontarget
                plot(i,100*R(i),'rs','markerface','r')
            end
        end
        xlim([0.5 length(freqs)+0.5])
        set(gca,'xtick',1:length(freqs),'xticklabel',freqs)
        xlabel('Frequency (kHz)')
        aa=axis;
        ylim([0 aa(4)])
    else
        contourf(str2num(Frqs), 60+str2num(Lvls),imgaussfilt(R,1.5))
        xlabel('Frequency (kHz)')
        ylabel('Level (dB SPL)')
        colormap(jet)
    end
    set(gca,'fontsize',8)
    title([{'Tone Responses'}])
    
else
    nospec=1;
    
end
set(gca,'box','off')

%% Bar plots of response rates
bhvStat.earlyCount = earlyCount;
bhvStat.hitCount = hitCount;
bhvStat.falseAlarmCount = falseAlarmCount;
bhvStat.totalTrials = totalTrials;

if phaseChoice  ~= 1 & ~nospec
    subplot(1,8,7)
elseif nospec
    subplot(1,8,6)
end
cla
bar(1,100*earlyCount/totalTrials,'facecolor','g','edgecolor','g');
hold on
bar(2,100*hitCount/totalTrials,'facecolor','b','edgecolor','none');
bar(3,100*falseAlarmCount/totalTrials,'facecolor','r','edgecolor','none');
title('Response Rates')
if phaseChoice ~=4 && ~isempty(levelVec)
    dprime = roundTo(norminv(hitCount/totalTrials) - norminv(earlyCount/totalTrials),2);
    set(gca,'xtick',1:2,'xticklabel',{'E','T'})
else
    dprime = roundTo(norminv(hitCount/totalTrials) - norminv(falseAlarmCount/totalTrials),2);
    set(gca,'xtick',1:3,'xticklabel',{'E','T','NT'})
end
bhvStat.dprime = dprime;
aa=axis;
if phaseChoice ~= 1
    text(0.25,aa(4)*1.1,['d''=' num2str(dprime)]);
end
set(gca,'fontsize',8)
ylim([0 aa(4)*1.2])
set(gca,'box','off')

%% Task Engagement histograms
Rtrials = find(sum(totalData,2)>0);
optTimesTemp = timeStamp(Rtrials,2);
formatIn = 'MM/dd/yyyy HH:mm:ss';
optTimes=datetime;
for i = 1:length(optTimesTemp)
    DateString = optTimesTemp{i};
    optTimes(i) = datetime(DateString,'InputFormat',formatIn);
end
for i = 1:length(optTimes)
    dates(i,:) = [month(optTimes(i)) day(optTimes(i))];
end
dates = unique(dates, 'Rows');
bhvStat.dates = dates';
bhvStat.phaseChoice = phaseChoice;
bhvStat.timeStamp = timeStamp;
bhvStat.totalData = totalData;
bhvStat.responseVec = responseVec;
if phaseChoice  ~= 1 & ~nospec
    subplot(2,8,16:16.5)
elseif nospec
    subplot(2,8,15:15.5)
end
cla
DoptTrial = round(days(optTimes-optTimes(1)))+1;
b=1:max(DoptTrial);
h=hist(DoptTrial,b,'facecolor','k','edgecolor','none');
bar(b,h,1,'facecolor','k','edgecolor','none');
if max(DoptTrial) == 1
    xlim([1 1.5])
else
    xlim([1 max(DoptTrial)])
end
bhvStat.DoptTrialsDays= [h; b];
B = unique(DoptTrial);
out = [B'  histc(DoptTrial,B)'];
bhvStat.DoptTrialsDaysData = out(:,2);
aa=axis;
xlabel('Training Day')
set(gca,'fontsize',8)
ylabel('Trials')
set(gca,'box','off')

if phaseChoice  ~= 1 & ~nospec
    subplot(2,8,8:8.5)
elseif nospec
    subplot(2,8,7:7.5)
end
cla
idxDays = round(days(optTimes-optTimes(1)))+1;
trialsPerMin=[];
for i = 1:length(b)
    try
        idx = find(idxDays == i);
        optTimesIdx=optTimes(idx);
        idxMin = minutes(optTimesIdx-optTimesIdx(1));
        tidx=[];
        for ii = 0:ceil(max(idxMin))-1
            tidx = [tidx; sum(idxMin>=ii & idxMin<ii+1)];
        end
        trialsPerMin = [trialsPerMin; mean(tidx)];
    catch
        trialsPerMin = [trialsPerMin; nan];
        
    end
end

plot(1:length(trialsPerMin),trialsPerMin)
if length(trialsPerMin) == 1
    xlim([1 1.5])
else
    xlim([1 length(trialsPerMin)])
end
xlabel('Training Day')
ylabel('Trials/Min.')
bhvStat.trialsPerMinData= trialsPerMin;
aa=axis;
title([{'Task Engagement'}])
set(gca,'fontsize',8)
set(gca,'box','off')

%% Creating standardized output structures
bhvOutput.earlyCount = earlyCount;
bhvOutput.falseAlarmCount = falseAlarmCount;
bhvOutput.hitCount = hitCount;
bhvOutput.levelVec = levelVec;
bhvOutput.totalData = totalData;
bhvOutput.phaseChoice = phaseChoice;
bhvOutput.responseVec = responseVec;
if ~exist('target', 'var')
    target = [];
end
bhvOutput.target = target;
if ~exist('nontarget', 'var')
    nontarget = [];
end
bhvOutput.nontarget = nontarget;
bhvOutput.timeStamp = timeStamp;
bhvOutput.toneVec = toneVec;
bhvOutput.totalTrials = totalTrials;