%runs monthly analysis on tonebox behavior and creates graphs to visualize
%the data
%IMPORTANT: cannot be run until "toneboxBehaviorSetup" script has been run
%to created compiled data files and monthly data files for the tonebox to
%be run through analysis
%Ilan Goldstein (11/2019)
numBoxes = input('Number of toneboxes to add to analysis: ');              %request user for number of boxes to be analyzed
for i = 1:numBoxes
    runBoxes(i) = input('Tonebox to be add to analysis: ');                %user input tonebox numbers for analysis
    boxNames{i} = strcat('tonebox',num2str(runBoxes(i)));
end
dataLoc = 'D:\Devices';                                                    %location of tonebox data, files should be separated into tonebox folders
genPath = 'C:\Users\PsiDev\Desktop\Ilan_Psignal';                          %location of analysis files (toneboxBehavior, toneBoxAnalysis, correctToneLevelFreq)
cd(dataLoc)
storedBoxes = dir('*tonebox*');
files = {};
%selecting files to add to analysis, different training phases can be included by uncommenting training phase for-loop%
for i = 1:numBoxes
    boxName = strcat('tonebox',num2str(runBoxes(i)));
    monthlyFiles = fullfile(dataLoc,boxName,'monthly data')
    cd(monthlyFiles)
    boxFiles = dir(monthlyFiles);
    count = 1;
%     for ii = 1:length(boxFiles)                                            %habituation
%         if contains(boxFiles(ii).name, 'habituation')                    
%             fileName = fullfile(dataLoc,boxName,'monthly data',boxFiles(ii).name);
%             files{count,i} = fileName;                                   
%             count = count + 1;                                           
%         end                                                              
%     end                                                                  
%     for ii = 1:length(boxFiles)                                            %shaping
%         if contains(boxFiles(ii).name, 'shaping')                        
%             fileName = fullfile(dataLoc,boxName,'monthly data',boxFiles(ii).name);
%             files{count,i} = fileName;                                   
%             count = count + 1;                                           
%         end                                                              
%     end                                                                  
    for ii = 1:length(boxFiles)                                            %detection
        if contains(boxFiles(ii).name, 'detection')                        
            fileName = fullfile(dataLoc,boxName,'monthly data',boxFiles(ii).name);
            files{count,i} = fileName;                                     
            count = count + 1;                                             
        end                                                                
    end                                                                    
%     for ii = 1:length(boxFiles)                                            %discrimination
%         if contains(boxFiles(ii).name, 'discrimination')                   
%             fileName = fullfile(dataLoc,boxName,'monthly data',boxFiles(ii).name);
%             files{count,i} = fileName;                                     
%             count = count + 1;                                             
%         end                                                                
%     end                                                                    
end
cd(genPath)
%running each monthly data file through "toneBoxAnalysisMonthly" to%
%generate behavior graphs and add monthly data to cross-monthly matrices%
for i = 1:size(files,2)  
    tpd = [];
    for ii = 1:size(files,1)
        if ~isempty(files{ii,i})
            loadFile = files{ii,i};
            [bhvOutput] = toneBoxAnalysisMonthly(loadFile);
            figTitle = extractAfter(loadFile,11);
            suptitle(figTitle)
            figName = extractAfter(figTitle,"\");
            figName = figName(1:end-4);
            figures = 'figures';
            set(gcf, 'WindowStyle', 'Docked')
            figSave = char(fullfile(dataLoc,boxNames(i),figures,figName));
            savefig(figSave)
            load(loadFile)
            start(ii,i) = bhvOutput.timeStamp(1);
            earlyRates(ii,i) = bhvOutput.earlyRate;
            hitRates(ii,i) = bhvOutput.hitRate;
            falseAlarmRates(ii,i) = bhvOutput.falseAlarmRate;
            dprimes(ii,i) = bhvOutput.dprime;
            targetResps{ii,i} = bhvOutput.targetResp;
            nontargetResps{ii,i} = bhvOutput.nontargetResp;
            t{ii,i} = bhvOutput.t;
            lickResponseE{ii,i} = bhvOutput.lickResponseE;
            lickResponseH{ii,i} = bhvOutput.lickResponseH;
            lickResponseF{ii,i} = bhvOutput.lickResponseF;
            tpd = [tpd; bhvOutput.out];
            fileNames{ii,i} = figName;
        end
    end
    trialsPerDay{i} = tpd;
end
tpdDate = {};
%calculating number of active trials per day for each monthly file%
for i = 1:size(trialsPerDay,2)
    new = [];
    new = find(trialsPerDay{i}(:,1) == 1);
    new(end+1) = size(trialsPerDay{i},1);
    newCount = 1;
    for ii = 1:size(trialsPerDay{i},1)
        if ii == new(newCount) && new(newCount) ~= new(end)
            tpdDate{i}(ii,1) = datetime(start{newCount,i},'InputFormat','MM/dd/yyyy HH:mm:ss','Format','MM-dd-yyyy');
            newCount = newCount + 1;
        else
            tpdDate{i}(ii,1) = datetime(start{newCount-1,i},'InputFormat','MM/dd/yyyy HH:mm:ss','Format','MM-dd-yyyy') + caldays(trialsPerDay{i}(ii,1)-1);
        end
    end
end    
%plotting monthly behavioral data onto single graphs for comparison across months%
for i = 1:size(files,2)           
    boxFigs(i) = figure;
    boxFigs(i).Name = strcat('tonebox_',num2str(i));
    colormap = jet;
    for ii = 1:size(files,1)
        if ~isempty(files{ii,i}) %&& ~contains(files{ii,i}, 'habituation') && ~contains(files{ii,i}, 'shaping')
            subplot(1,8,1)
            title('trial licking')
            hold on
            plot(targetResps{ii,i}(:,1),targetResps{ii,i}(:,2),'DisplayName',char(fileNames{ii,i}),'Color',colormap(ii*5,:))
            hold off
            respRates(ii,:,i) = [earlyRates(ii,i) hitRates(ii,i) falseAlarmRates(ii,i)];
            ylabel('percent (%)')
            xlabel('time (s)')
            subplot(1,8,2)
            title({'early', 'response latency'})
            LickSum = (lickResponseE{ii,i} + lickResponseH{ii,i} + lickResponseF{ii,i});
            ELR = 100*lickResponseE{ii,i}/sum(LickSum);
            sELR = smooth(ELR,10);
            Eidx = find(t{ii,i}<=1);
            hold on
            area(t{ii,i}(Eidx),sELR(Eidx),'facecolor',colormap(ii*5,:),'facealpha',.3,'edgecolor','none','DisplayName',char(fileNames{ii,i}));
            hold off
            xlabel('time (s)')
            subplot(1,8,3)
            title({'hit', 'response latency'})
            HLR = 100*lickResponseH{ii,i}/sum(LickSum);
            sHLR = smooth(HLR,10);
            idx = find(t{ii,i}>=1);
            idx = [min(idx)-1 idx];
            hold on
            area(t{ii,i}(idx),sHLR(idx),'facecolor',colormap(ii*5,:),'facealpha',.3,'edgecolor','none','DisplayName',char(fileNames{ii,i}));
            hold off
            xlabel('time (s)')
            subplot(1,8,4)
            title({'false alarm', 'response latency'})
            FLR = 100*lickResponseF{ii,i}/sum(LickSum);
            sFLR = smooth(FLR,10);
            hold on
            area(t{ii,i}(idx),sFLR(idx),'facecolor',colormap(ii*5,:),'facealpha',.3,'edgecolor','none','DisplayName',char(fileNames{ii,i}));
            hold off
            xlabel('time (s)')
            dprLabel(ii) = dprimes(ii,i);
        end
    end
    leg1 = legend();
    set(leg1,'Interpreter','none')
    respBar = respRates(any(respRates(:,:,i),2),:,i);
    subplot(1,8,5:6)
    title('response rates')
    hold on
    b = bar(respBar);
    xtips = 1:length(dprLabel);
    ytips = respBar(:,2);
    labels = num2str(dprLabel);
    text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
    leg2 = legend('early', 'target', 'nontarget');
    xticks([1:length(respBar)])
    xticklabels(leg1.String)
    set(gca,'TickLabelInterpreter','none')
    xtickangle(-15)
    titleName = strcat(boxNames(i),[' monthly average behavior']);
    suptitle(titleName)
    clear respBar dprLabel
    subplot(1,8,7:8)
    bar(tpdDate{i},trialsPerDay{i}(:,2),1)
    title('task engagement')
    ylabel('trials per day')
    set(gcf, 'WindowStyle', 'Docked')
    figName = 'month_avg_behavior_performance';
    figSave = char(fullfile(dataLoc,boxNames(i),figures,figName));
    savefig(figSave)
end
