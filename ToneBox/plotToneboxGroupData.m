function GroupStats = plotToneboxGroupData(phase, lickOGramTarget, lickOGramNonTarget, earlyCount, hitCount, falseAlarmCount, ...
    totalTrials, LatencyH, LatencyF, LatencyE, DoptTrialsDays, D, DoptTrialsDaysCages,DoptTrialsPerMinute)
%plotToneboxGroupData polots group data for Francis et al. 2019
%Nikolas Francis, 12/2018

fignum= 10000*round(rand,4)
figure(fignum)
set(gcf,'position',[14 649 1884 311])

%Time limits per trial
if phase == 1
    Tidx = [0 10];
else
    Tidx = [0 4];
end

%% Plot the lick-o-gram, ie., the average lick response
figure(fignum)
subplot(1,4,1)
cla

if phase ==4
    plot(0,0,'b')
    hold on
    plot(0,0,'r')
    plot(0,0,'g')
    legend('T','NT','E','autoupdate','off')
elseif phase ~= 1
    plot(0,0,'b')
    hold on
    plot(0,0,'g')
    legend('T','E','autoupdate','off')
else
    plot(0,0,'b')
    hold on
    legend('T','autoupdate','off')
end

mulickOGramTarget = nanmean(lickOGramTarget,2);
stelickOGramTarget = 2*std(lickOGramTarget,[],2)./sqrt(size(lickOGramTarget,2));
shadedErrorBar(linspace(0,Tidx(2),length(mulickOGramTarget)),mulickOGramTarget,stelickOGramTarget,'b',1);
if ~isempty(lickOGramNonTarget)
    mulickOGramNonTarget = nanmean(lickOGramNonTarget,2);
    stelickOGramNonTarget = 2*std(lickOGramNonTarget,[],2)./sqrt(size(lickOGramNonTarget,2));
    shadedErrorBar(linspace(0,Tidx(2),length(mulickOGramNonTarget)),mulickOGramNonTarget,stelickOGramNonTarget,'b',1);
end
legend boxoff
xlim([0 Tidx(2)]);
xlabel('Time(s)')
hold on
title('Lick-o-gram')
set(gca,'fontsize',8)
aa=axis;
ylim([aa(3) aa(4)])
if phase ~=1
    %Box around target
    y1 = [0 aa(4)]*.5;
    x1 = [1 1];
    y2 = [aa(4) aa(4)]*.5;
    x2 = [2 2];
    area([x1 x2],[y1 y2],'facecolor','k','facealpha',.25,'edgecolor','none');
end
ylabel('Percent (%)')
set(gca,'box','off')

%% Plot response latency, aka first lick
figure(fignum)
subplot(1,4,2)
cla

% LatencyH = [zeros(size(lickOGramTarget,1)-size(LatencyH,1),size(LatencyH,2)); LatencyH];
% muLatencyH = nanmean(LatencyH,2);
% steLatencyH= 2*std(LatencyH,[],2)./sqrt(size(LatencyH,2));
% shadedErrorBar(linspace(0,Tidx(2),length(muLatencyH)),muLatencyH,steLatencyH,'b',1);
hold on
if ~isempty(LatencyF)
    LatencyF = [zeros(size(lickOGramNonTarget,1)-size(LatencyF,1),size(LatencyF,2)); LatencyF];
    muLatencyF = nanmean(LatencyF,2);
    steLatencyF = 2*std(LatencyF,[],2)./sqrt(size(LatencyF,2));
    shadedErrorBar(linspace(0,Tidx(2),length(muLatencyF)),muLatencyF,steLatencyF,'r',1);
    
        Latency=[];
    for i = 1:size(LatencyE,2)
        Latency = [Latency smooth([LatencyE(:,i); LatencyF(:,i)],10)];
    end
    muLatencyE = nanmean(Latency(1:size(LatencyE,1),:),2);
    steLatencyE = 2*std(Latency(1:size(LatencyE,1),:),[],2)./sqrt(size(Latency,2));
    shadedErrorBar(linspace(0,1,length(muLatencyE)),muLatencyE,steLatencyE,'g',1);
    hold on
    muLatencyF = nanmean(Latency(size(LatencyE,1)+1:size(LatencyE,1)+size(LatencyF,1),:),2);
    steLatencyF = 2*std(Latency(size(LatencyE,1)+1:size(LatencyE,1)+size(LatencyF,1),:),[],2)./sqrt(size(Latency,2));
    shadedErrorBar(linspace(1,4,length(muLatencyF)),muLatencyF,steLatencyF,'b',1);

end
if ~isempty(LatencyE)
    Latency=[];
    for i = 1:size(LatencyE,2)
        Latency = [Latency smooth([LatencyE(:,i); LatencyH(:,i)],10)];
    end
    muLatencyE = nanmean(Latency(1:size(LatencyE,1),:),2);
    steLatencyE = 2*std(Latency(1:size(LatencyE,1),:),[],2)./sqrt(size(Latency,2));
    shadedErrorBar(linspace(0,1,length(muLatencyE)),muLatencyE,steLatencyE,'g',1);
    hold on
    muLatencyH = nanmean(Latency(size(LatencyE,1)+1:size(LatencyE,1)+size(LatencyH,1),:),2);
    steLatencyH = 2*std(Latency(size(LatencyE,1)+1:size(LatencyE,1)+size(LatencyH,1),:),[],2)./sqrt(size(Latency,2));
    shadedErrorBar(linspace(1,4,length(muLatencyH)),muLatencyH,steLatencyH,'b',1);
end

xlim([0 Tidx(2)]);
hold on
title('Response Latency')
set(gca,'fontsize',8)
aa=axis;
ylim([aa(3) aa(4)])
if phase ~=1
    %Box around target
    y1 = [0 aa(4)]*.5;
    x1 = [1 1];
    y2 = [aa(4) aa(4)]*.5;
    x2 = [2 2];
    area([x1 x2],[y1 y2],'facecolor','k','facealpha',.25,'edgecolor','none');
end
set(gca,'box','off')

%Stats
t=linspace(0,Tidx(2),size(Latency,1));
[j m] = max(Latency);
GroupStats.Latency.peak = t(m);
GroupStats.Latency.peakmu = mean(t(m));
GroupStats.Latency.peakste= std(t(m))./sqrt(length(m));
[GroupStats.Latency.h GroupStats.Latency.p] = ttest(t(m),1);

%Plot response accuracy via latency
figure(round(fignum/2))
subplot(1,2,1)
cla

muH=[];
muF=[];
muE=[];
muH= 100*(sum(Latency(20:end,:)./sum(Latency)));
if ~isempty(falseAlarmCount)
    muF= 100*(sum(Latency(20:end,:)./sum(Latency)));
end
if ~isempty(earlyCount)
    muE= 100*(sum(Latency(1:20,:)./sum(Latency)));
end
if phase == 1
    boxplot(muH, 'Notch','off','Labels', {'W'},'Whisker',1.5)
elseif phase == 2 || phase == 3
    boxplot([muH' muE'], 'Notch','off','Labels', {'T','E'},'Whisker',1.5)
else
    boxplot([muH' muF' muE'], 'Notch','off','Labels', {'T','F','E'},'Whisker',1.5)
end
hold on
aa=axis;
if ~isempty(falseAlarmCount)
    dprime = roundTo(norminv(muH./100) - norminv(muF./100),2);
    text(0.25,aa(4)*1.1,['d''=' num2str(dprime)]);
end
if ~isempty(earlyCount)
    dprime = mean(roundTo(norminv(muH./100) - norminv(muE./100),2));
    text(aa(1),aa(4)*1.1,['d''=' num2str(dprime)]);
end
aa=axis;
set(gca,'fontsize',8)
ylim([0 aa(4)*1.2])
set(gca,'box','off')
title('Accuracy: Response Latency')

%Stats
GroupStats.Latency.H = muH;
GroupStats.Latency.Hmu = mean(muH);
GroupStats.Latency.Hste= std(muH)./sqrt(length(muH));
GroupStats.Latency.F= muF;
GroupStats.Latency.Fmu = mean(muF);
GroupStats.Latency.Fste= std(muF)./sqrt(length(muF));
GroupStats.Latency.E= muE;
GroupStats.Latency.Emu = mean(muE);
GroupStats.Latency.Este= std(muE)./sqrt(length(muE));

%% Plot response rates
figure(round(fignum./2))
subplot(1,2,2)
cla

muH=[];
muF=[];
muE=[];
muH= 100*(hitCount./totalTrials);
if ~isempty(falseAlarmCount)
    muF = 100*(falseAlarmCount./totalTrials);
end
if ~isempty(earlyCount)
    muE = 100*(earlyCount./totalTrials);
end
if phase == 1
    boxplot(muH, 'Notch','off','Labels', {'W'},'Whisker',1.5)
elseif phase == 2 || phase == 3
    boxplot([muH' muE'], 'Notch','off','Labels', {'T','E'},'Whisker',1.5)
else
    boxplot([muH' muF' muE'], 'Notch','off','Labels', {'T','F','E'},'Whisker',1.5)
end
hold on
aa=axis;
if ~isempty(falseAlarmCount)
    dprime = roundTo(norminv(muH./100) - norminv(muF./100),2);
    text(0.25,aa(4)*1.1,['d''=' num2str(dprime)]);
end
if ~isempty(earlyCount)
    dprime = mean(roundTo(norminv(muH./100) - norminv(muE./100),2));
    text(aa(1),aa(4)*1.1,['d''=' num2str(dprime)]);
end
aa=axis;
set(gca,'fontsize',8)
ylim([0 aa(4)*1.2])
set(gca,'box','off')
title('Responsiveness')

%Stats
GroupStats.ResponseRates.H = muH;
GroupStats.ResponseRates.E = muE;
GroupStats.ResponseRates.F = muF;

%% Task engagement
figure(fignum)
subplot(1,4,3)
cla

muDoptTrialsPerMinute = nanmean(DoptTrialsPerMinute,2);
steDoptTrialsPerMinute= 2*nanstd(DoptTrialsPerMinute,[],2)./sqrt(size(DoptTrialsPerMinute,2));
shadedErrorBar(D,muDoptTrialsPerMinute,steDoptTrialsPerMinute,'k',1);
xlim([1 size(DoptTrialsDays,1)])
ylabel('Trials/Min.')
title([{'Task Engagement'}])
set(gca,'fontsize',8)
set(gca,'box','off')

subplot(1,4,4)
cla

muDoptTrialsDays = nanmean(DoptTrialsDays,2);
steDoptTrialsDays= 2*std(DoptTrialsDays,[],2)./sqrt(size(DoptTrialsDays,2));;
shadedErrorBar(D,muDoptTrialsDays,steDoptTrialsDays,'k',1);
xlabel('Training Day')
set(gca,'fontsize',8)
ylabel('Trials')
xlim([1 size(DoptTrialsDays,1)])
set(gca,'box','off')

%Stats
GroupStats.TaskEng.DoptTrialsDays = DoptTrialsDaysCages;
GroupStats.TaskEng.DoptTrialsPerMinute = DoptTrialsPerMinute;
GroupStats.TaskEng.muDoptTrialsDays = mean(DoptTrialsDaysCages);
GroupStats.TaskEng.muDoptTrialsPerMinute = muDoptTrialsPerMinute
GroupStats.TaskEng.steDoptTrialsDays = std(DoptTrialsDaysCages)./sqrt(length(DoptTrialsDaysCages));
GroupStats.TaskEng.steDoptTrialsPerMinute = steDoptTrialsPerMinute./2;
% 


