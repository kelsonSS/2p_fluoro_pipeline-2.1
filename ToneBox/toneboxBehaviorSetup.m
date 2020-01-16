% toneboxBehavior separates tonebox behavior data by tonebox, training
% phase, and month, then combines these data into single files to be
% analyzed using toneBoxAnalysis (Francis, 2019)
% Ilan Goldstein 08/2019
numBoxes = input('Number of toneboxes to add to analysis: ');              %request user for number of boxes to be analyzed
for i = 1:numBoxes
    runBoxes(i) = input('Tonebox to be add to analysis: ');                %user input tonebox numbers for analysis
    boxNames{i} = strcat('tonebox',num2str(runBoxes(i)));
end
dataLoc = 'D:\Devices';                                                    %location of tonebox data, files should be separated into tonebox folders
genPath = 'C:\Users\PsiDev\Desktop\Ilan_Psignal';                          %location of analysis files (toneboxBehaviorSetup, toneBoxAnalysis, correctToneLevelFreq, toneboxBehaviorMonthly)
cd(dataLoc)
storedBoxes = dir('*tonebox*');
files = {};
compData = {};
%accessing behavior data files and creating structure arrays for all files
%corrsponding to each tonebox with usable data that are grouped together 
%into the cell array "compData"
tic
for i = 1:numBoxes
    boxName = strcat('tonebox',num2str(runBoxes(i)));
    cd(dataLoc)
    boxFiles = dir(boxName);
    count = 1;
    for ii = 1:length(boxFiles)
        if contains(boxFiles(ii).name, 'performance') && ~contains(boxFiles(ii).name, 'stop_error') && ~contains(boxFiles(ii).name, 'backup')
            loadFile = fullfile(dataLoc,boxName,boxFiles(ii).name);
            files = [files; loadFile];
            cd(genPath)
            [bhvStat bhvOutput] = toneBoxAnalysis(loadFile);
            if ~isempty(bhvOutput)
                bhvData(count) = bhvOutput;
                count = count + 1;
            end
            set(gcf, 'WindowStyle', 'Docked')
        else
        end
    end
    compData{i} = bhvData;
    clear bhvData bhvOutput
end
%find start time for each data file within each tonebox cell
for i = 1:length(compData)
    for ii = 1:length(compData{i})
        start{i}(ii,1) = [datetime(compData{i}(ii).timeStamp{1,2})];
    end
end
%order start times of files within each tonebox cell
for i = 1:length(start)
    dateOrder{i} = sort(start{i});
end
%index initial file order and match to ordered file dates
for i = 1:length(compData)
    for ii = 1:length(dateOrder{1,i})
        ordIDX = find(start{i} == dateOrder{1,i}(ii));
        dateOrder{2,i}(ii) = ordIDX(1);
    end
end
%reorder files in each tonebox cell according to file start date order
cmpOrdrData = {};
for i = 1:length(compData)
    for ii = 1:length(dateOrder{2,i})
        cmpOrdrData{i}(ii) = compData{i}(dateOrder{2,i}(ii));
    end
end
%separate files by tonebox based on training phase and create cell arrays for each phase 
shapeData = {};
habitData = {};
detectData = {};
discrimData = {};
for i = 1:length(cmpOrdrData)
    shapeCount = 1;
    habitCount = 1;
    detectCount = 1;
    discrimCount = 1;
    for ii = 1:length(cmpOrdrData{i})
        if cmpOrdrData{i}(ii).totalTrials >= 1000                          %only data files with >=1000 trials are kept at this stage
            if cmpOrdrData{i}(ii).phaseChoice == 1
                habitData{i}(habitCount) = cmpOrdrData{i}(ii);
                habitCount = habitCount + 1;
            elseif cmpOrdrData{i}(ii).phaseChoice == 2
                shapeData{i}(shapeCount) = cmpOrdrData{i}(ii);
                shapeCount = shapeCount + 1;
            elseif cmpOrdrData{i}(ii).phaseChoice == 3
                detectData{i}(detectCount) = cmpOrdrData{i}(ii);
                detectCount = detectCount + 1;
            elseif cmpOrdrData{i}(ii).phaseChoice == 4
                discrimData{i}(discrimCount) = cmpOrdrData{i}(ii);
                discrimCount = discrimCount + 1;
            end
        else
        end
    end
end
toc
%combining habituation behavior data from different files into one structure for each tonebox
habitComp = {};
for i = 1:length(habitData)
    earlyCount = 0;
    falseAlarmCount = 0;
    hitCount = 0;
    levelVec = [];
    totalData = [];
    phaseChoice = [];
    responseVec = [];
    target = [];
    nontarget = [];
    timeStamp = [];
    toneVec = [];
    totalTrials = 0;
    if ~isempty(habitData{i})
        for ii = 1:length(habitData{i})
            earlyCount = earlyCount + habitData{i}(ii).earlyCount;
            falseAlarmCount = falseAlarmCount + habitData{i}(ii).falseAlarmCount;
            hitCount = hitCount + habitData{i}(ii).hitCount;
            levelVec = [levelVec habitData{i}(ii).levelVec];
            totalData = [totalData; habitData{i}(ii).totalData];
            phaseChoice = habitData{i}(ii).phaseChoice;
            responseVec = [responseVec habitData{i}(ii).responseVec];
            target = habitData{i}(ii).target;
            nontarget = habitData{i}(ii).nontarget;
            timeStamp = [timeStamp; habitData{i}(ii).timeStamp];
            toneVec = [toneVec habitData{i}(ii).toneVec];
            totalTrials = totalTrials + habitData{i}(ii).totalTrials;
        end
        habitComp{i}.earlyCount = earlyCount;
        habitComp{i}.falseAlarmCount = falseAlarmCount;
        habitComp{i}.hitCount = hitCount;
        habitComp{i}.levelVec = levelVec;
        habitComp{i}.totalData = totalData;
        habitComp{i}.phaseChoice = phaseChoice;
        habitComp{i}.responseVec = responseVec;
        habitComp{i}.target = target;
        habitComp{i}.nontarget = nontarget;
        habitComp{i}.timeStamp = timeStamp;
        habitComp{i}.toneVec = toneVec;
        habitComp{i}.totalTrials = totalTrials;
    else
    end
end
%combining shaping behavior data from different files into one structure for each tonebox
shapeComp = {};
for i = 1:length(shapeData)
    earlyCount = 0;
    falseAlarmCount = 0;
    hitCount = 0;
    levelVec = [];
    totalData = [];
    phaseChoice = [];
    responseVec = [];
    target = [];
    nontarget = [];
    timeStamp = [];
    toneVec = [];
    totalTrials = 0;
    if ~isempty(shapeData{i})
        for ii = 1:length(shapeData{i})
            earlyCount = earlyCount + shapeData{i}(ii).earlyCount;
            falseAlarmCount = falseAlarmCount + shapeData{i}(ii).falseAlarmCount;
            hitCount = hitCount + shapeData{i}(ii).hitCount;
            levelVec = [levelVec shapeData{i}(ii).levelVec];
            totalData = [totalData; shapeData{i}(ii).totalData];
            phaseChoice = shapeData{i}(ii).phaseChoice;
            responseVec = [responseVec shapeData{i}(ii).responseVec];
            target = shapeData{i}(ii).target;
            nontarget = shapeData{i}(ii).nontarget;
            timeStamp = [timeStamp; shapeData{i}(ii).timeStamp];
            toneVec = [toneVec shapeData{i}(ii).toneVec];
            totalTrials = totalTrials + shapeData{i}(ii).totalTrials;
        end
        shapeComp{i}.earlyCount = earlyCount;
        shapeComp{i}.falseAlarmCount = falseAlarmCount;
        shapeComp{i}.hitCount = hitCount;
        shapeComp{i}.levelVec = levelVec;
        shapeComp{i}.totalData = totalData;
        shapeComp{i}.phaseChoice = phaseChoice;
        shapeComp{i}.responseVec = responseVec;
        shapeComp{i}.target = target;
        shapeComp{i}.nontarget = nontarget;
        shapeComp{i}.timeStamp = timeStamp;
        shapeComp{i}.toneVec = toneVec;
        shapeComp{i}.totalTrials = totalTrials;
    else
    end
end
%combining detection behavior data from different files into one structure for each tonebox
detectComp = {};
for i = 1:length(detectData)
    earlyCount = 0;
    falseAlarmCount = 0;
    hitCount = 0;
    levelVec = [];
    totalData = [];
    phaseChoice = [];
    responseVec = [];
    target = [];
    nontarget = [];
    timeStamp = [];
    toneVec = [];
    totalTrials = 0;
    formatIDX = [];
    if ~isempty(detectData{i})
        formats = input('How many data formats exist? ')                   %used for detection files with multiple stimulus types
        for ii = 1:(formats)
            formatIDX(ii) = input('File index after end of format: ')      %found under "detectData", 1st = field # starting 2nd stimulus type, last = field # after last data containing field
        end
        formatNum = 1;
        for ii = 1:length(detectData{i})
            if ii < formatIDX(formatNum)
                earlyCount = earlyCount + detectData{i}(ii).earlyCount;
                falseAlarmCount = falseAlarmCount + detectData{i}(ii).falseAlarmCount;
                hitCount = hitCount + detectData{i}(ii).hitCount;
                levelVec = [levelVec detectData{i}(ii).levelVec];
                totalData = [totalData; detectData{i}(ii).totalData];
                phaseChoice = detectData{i}(ii).phaseChoice;
                responseVec = [responseVec detectData{i}(ii).responseVec];
                target = detectData{i}(ii).target;
                nontarget = detectData{i}(ii).nontarget;
                timeStamp = [timeStamp; detectData{i}(ii).timeStamp];
                toneVec = [toneVec detectData{i}(ii).toneVec];
                totalTrials = totalTrials + detectData{i}(ii).totalTrials;
            elseif ii == formatIDX(formatNum)
                earlyCount = earlyCount + detectData{i}(ii).earlyCount;
                falseAlarmCount = falseAlarmCount + detectData{i}(ii).falseAlarmCount;
                hitCount = hitCount + detectData{i}(ii).hitCount;
                levelVec = [levelVec detectData{i}(ii).levelVec];
                totalData = [totalData; detectData{i}(ii).totalData];
                phaseChoice = detectData{i}(ii).phaseChoice;
                responseVec = [responseVec detectData{i}(ii).responseVec];
                target = detectData{i}(ii).target;
                nontarget = detectData{i}(ii).nontarget;
                timeStamp = [timeStamp; detectData{i}(ii).timeStamp];
                toneVec = [toneVec detectData{i}(ii).toneVec];
                totalTrials = totalTrials + detectData{i}(ii).totalTrials;
                detectComp{i}(formatNum).earlyCount = earlyCount;
                detectComp{i}(formatNum).falseAlarmCount = falseAlarmCount;
                detectComp{i}(formatNum).hitCount = hitCount;
                detectComp{i}(formatNum).levelVec = levelVec;
                detectComp{i}(formatNum).totalData = totalData;
                detectComp{i}(formatNum).phaseChoice = phaseChoice;
                detectComp{i}(formatNum).responseVec = responseVec;
                detectComp{i}(formatNum).target = target;
                detectComp{i}(formatNum).nontarget = nontarget;
                detectComp{i}(formatNum).timeStamp = timeStamp;
                detectComp{i}(formatNum).toneVec = toneVec;
                detectComp{i}(formatNum).totalTrials = totalTrials;
                formatNum = formatNum + 1;
                earlyCount = 0;
                falseAlarmCount = 0;
                hitCount = 0;
                levelVec = [];
                totalData = [];
                phaseChoice = [];
                responseVec = [];
                target = [];
                nontarget = [];
                timeStamp = [];
                toneVec = [];
                totalTrials = 0;
            end
        end
    else
    end
end
%combining discrimination behavior data from different files into one structure for each tonebox
discrimComp = {};
for i = 1:length(discrimData)
    earlyCount = 0;
    falseAlarmCount = 0;
    hitCount = 0;
    levelVec = [];
    totalData = [];
    phaseChoice = [];
    responseVec = [];
    target = [];
    nontarget = [];
    timeStamp = [];
    toneVec = [];
    totalTrials = 0;
    if ~isempty(discrimData{i})
        for ii = 1:length(discrimData{i})
            earlyCount = earlyCount + discrimData{i}(ii).earlyCount;
            falseAlarmCount = falseAlarmCount + discrimData{i}(ii).falseAlarmCount;
            hitCount = hitCount + discrimData{i}(ii).hitCount;
            levelVec = [levelVec discrimData{i}(ii).levelVec];
            totalData = [totalData; discrimData{i}(ii).totalData];
            phaseChoice = discrimData{i}(ii).phaseChoice;
            responseVec = [responseVec discrimData{i}(ii).responseVec];
            target = discrimData{i}(ii).target;
            nontarget = discrimData{i}(ii).nontarget;
            timeStamp = [timeStamp; discrimData{i}(ii).timeStamp];
            toneVec = [toneVec discrimData{i}(ii).toneVec];
            totalTrials = totalTrials + discrimData{i}(ii).totalTrials;
        end
        discrimComp{i}.earlyCount = earlyCount;
        discrimComp{i}.falseAlarmCount = falseAlarmCount;
        discrimComp{i}.hitCount = hitCount;
        discrimComp{i}.levelVec = levelVec;
        discrimComp{i}.totalData = totalData;
        discrimComp{i}.phaseChoice = phaseChoice;
        discrimComp{i}.responseVec = responseVec;
        discrimComp{i}.target = target;
        discrimComp{i}.nontarget = nontarget;
        discrimComp{i}.timeStamp = timeStamp;
        discrimComp{i}.toneVec = toneVec;
        discrimComp{i}.totalTrials = totalTrials;
    else
    end
end
%save combined phase data into single files for each tonebox
compDataLoc = 'compiled Data';
for i = 1:numBoxes
    habitName = 'habituationData';
    habitSave = fullfile(dataLoc,boxNames{i},compDataLoc,habitName);
    habitStruct = habitComp{i};
    save(habitSave,'habitStruct')
    shapeName = 'ShapingData';
    shapeSave = fullfile(dataLoc,boxNames{i},compDataLoc,shapeName);
    shapeStruct = shapeComp{i};
    save(shapeSave,'shapeStruct')
    detectName = 'detectionData';
    detectSave = fullfile(dataLoc,boxNames{i},compDataLoc,detectName);
    detectStruct = detectComp{i};
    save(detectSave,'detectStruct')
    discrimName = 'discriminationData';
    discrimSave = fullfile(dataLoc,boxNames{i},compDataLoc,discrimName);
    discrimStruct = discrimComp{i};
    save(discrimSave,'discrimStruct')
end
%separate combined habituation data by month into separate structures for each tonebox
tic
for i = 1:length(habitComp)
    if ~isempty(habitComp{i})
        for ii = 1:length(habitComp{i}.timeStamp)
            habitDate(ii,1) = month(habitComp{i}.timeStamp(ii,2));
        end
        habitUni = unique(habitDate,'stable');
        for ii = 1:length(habitUni)
            count = 1;
            for iii = 1:length(habitComp{i}.timeStamp)
                if month(habitComp{i}.timeStamp(iii,2)) == habitUni(ii)
                    habitMonth{i}(ii).levelVec(count) = habitComp{i}.levelVec(iii);
                    habitMonth{i}(ii).totalData(count,:) = habitComp{i}.totalData(iii,:);
                    habitMonth{i}(ii).responseVec(count) = habitComp{i}.responseVec(iii);
                    habitMonth{i}(ii).timeStamp(count,1) = habitComp{i}.timeStamp(iii,2);
                    habitMonth{i}(ii).toneVec(count) = habitComp{i}.toneVec(iii);
                    count = count + 1;
                else
                end
            end
            habitMonth{i}(ii).phaseChoice = habitComp{i}.phaseChoice;
            habitMonth{i}(ii).target = habitComp{i}.target;
            habitMonth{i}(ii).nontarget = habitComp{i}.nontarget;
        end
        clear habitDate habitUni
    end
end
%separate combined shaping data by month into separate structures for each tonebox
for i = 1:length(shapeComp)
    if ~isempty(shapeComp{i})
        for ii = 1:length(shapeComp{i}.timeStamp)
            shapeDate(ii,1) = month(shapeComp{i}.timeStamp(ii,2));
        end
        shapeUni = unique(shapeDate,'stable');
        for ii = 1:length(shapeUni)
            count = 1;
            for iii = 1:length(shapeComp{i}.timeStamp)
                if month(shapeComp{i}.timeStamp(iii,2)) == shapeUni(ii)
                    shapeMonth{i}(ii).levelVec(count) = shapeComp{i}.levelVec(iii);
                    shapeMonth{i}(ii).totalData(count,:) = shapeComp{i}.totalData(iii,:);
                    shapeMonth{i}(ii).responseVec(count) = shapeComp{i}.responseVec(iii);
                    shapeMonth{i}(ii).timeStamp(count,1) = shapeComp{i}.timeStamp(iii,2);
                    shapeMonth{i}(ii).toneVec(count) = shapeComp{i}.toneVec(iii);
                    count = count + 1;
                else
                end
            end
            shapeMonth{i}(ii).phaseChoice = shapeComp{i}.phaseChoice;
            shapeMonth{i}(ii).target = shapeComp{i}.target;
            shapeMonth{i}(ii).nontarget = shapeComp{i}.nontarget;
        end
        clear shapeDate shapeUni
    end
end
%separate combined detection data by month into separate structures for each tonebox
for i = 1:length(detectComp)
    if ~isempty(detectComp{i})
        for ii = 1:length(detectComp{i})
            for iii = 1:length(detectComp{i}(ii).timeStamp)
                detectDate(iii,1) = month(detectComp{i}(ii).timeStamp(iii,2));
            end
            detectUni = unique(detectDate,'stable');
            for iii = 1:length(detectUni)
                count = 1;
                for iv = 1:length(detectComp{i}(ii).timeStamp)
                    if month(detectComp{i}(ii).timeStamp(iv,2)) == detectUni(iii)
                        detectMonth{ii,i}(iii).levelVec(count) = detectComp{i}(ii).levelVec(iv);
                        detectMonth{ii,i}(iii).totalData(count,:) = detectComp{i}(ii).totalData(iv,:);
                        detectMonth{ii,i}(iii).responseVec(count) = detectComp{i}(ii).responseVec(iv);
                        detectMonth{ii,i}(iii).timeStamp(count,1) = detectComp{i}(ii).timeStamp(iv,2);
                        detectMonth{ii,i}(iii).toneVec(count) = detectComp{i}(ii).toneVec(iv);
                        count = count + 1;
                    else
                    end
                end
                detectMonth{ii,i}(iii).phaseChoice = detectComp{i}(ii).phaseChoice;
                detectMonth{ii,i}(iii).target = detectComp{i}(ii).target;
                detectMonth{ii,i}(iii).nontarget = detectComp{i}(ii).nontarget;
            end
        end
        clear detectDate detectUni
    end
end
%separate combined discrimination data by month into separate structures for each tonebox
for i = 1:length(discrimComp)
    if ~isempty(discrimComp{i})
        for ii = 1:length(discrimComp{i}.timeStamp)
            discrimDate(ii,1) = month(discrimComp{i}.timeStamp(ii,2));
        end
        discrimUni = unique(discrimDate,'stable');
        for ii = 1:length(discrimUni)
            count = 1;
            for iii = 1:length(discrimComp{i}.timeStamp)
                if month(discrimComp{i}.timeStamp(iii,2)) == discrimUni(ii)
                    discrimMonth{i}(ii).levelVec(count) = discrimComp{i}.levelVec(iii);
                    discrimMonth{i}(ii).totalData(count,:) = discrimComp{i}.totalData(iii,:);
                    discrimMonth{i}(ii).responseVec(count) = discrimComp{i}.responseVec(iii);
                    discrimMonth{i}(ii).timeStamp(count,1) = discrimComp{i}.timeStamp(iii,2);
                    discrimMonth{i}(ii).toneVec(count) = discrimComp{i}.toneVec(iii);
                    count = count + 1;
                else
                end
            end
            discrimMonth{i}(ii).phaseChoice = discrimComp{i}.phaseChoice;
            discrimMonth{i}(ii).target = discrimComp{i}.target;
            discrimMonth{i}(ii).nontarget = discrimComp{i}.nontarget;
        end
        clear discrimDate discrimUni
    end
end
toc
%save monthly habituation data into separate files by tonebox
clear levelVec totalData responseVec timeStamp toneVec phaseChoice target nontarget
for i = 1:length(habitMonth)
    if ~isempty(habitMonth{i})
        for ii = 1:length(habitMonth{i})
            date = [num2str(month(habitMonth{i}(ii).timeStamp(1))), '_', num2str(year(habitMonth{i}(ii).timeStamp(1)))];
            name = strcat(date,'_habituation');
            fileName = fullfile(dataLoc,boxNames{i},compDataLoc,name);
            levelVec = habitMonth{i}(ii).levelVec;
            totalData = habitMonth{i}(ii).totalData;
            responseVec = habitMonth{i}(ii).responseVec;
            timeStamp = habitMonth{i}(ii).timeStamp;
            toneVec = habitMonth{i}(ii).toneVec;
            phaseChoice = habitMonth{i}(ii).phaseChoice;
            target = habitMonth{i}(ii).target;
            nontarget = habitMonth{i}(ii).nontarget;
            save(fileName,'levelVec','totalData','responseVec','timeStamp','toneVec','phaseChoice','target','nontarget')
            clear date name fileName levelVec totalData responseVec timeStamp toneVec phaseChoice target nontarget
        end
    end
end
%save monthly shaping data into separate files by tonebox
for i = 1:length(shapeMonth)
    if ~isempty(shapeMonth{i})
        for ii = 1:length(shapeMonth{i})
            date = [num2str(month(shapeMonth{i}(ii).timeStamp(1))), '_', num2str(year(shapeMonth{i}(ii).timeStamp(1)))];
            name = strcat(date,'_shaping');
            fileName = fullfile(dataLoc,boxNames{i},compDataLoc,name);
            levelVec = shapeMonth{i}(ii).levelVec;
            totalData = shapeMonth{i}(ii).totalData;
            responseVec = shapeMonth{i}(ii).responseVec;
            timeStamp = shapeMonth{i}(ii).timeStamp;
            toneVec = shapeMonth{i}(ii).toneVec;
            phaseChoice = shapeMonth{i}(ii).phaseChoice;
            target = shapeMonth{i}(ii).target;
            nontarget = shapeMonth{i}(ii).nontarget;
            save(fileName,'levelVec','totalData','responseVec','timeStamp','toneVec','phaseChoice','target','nontarget')
            clear date name fileName levelVec totalData responseVec timeStamp toneVec phaseChoice target nontarget
        end
    end
end
%save monthly detection data into separate files by tonebox
for i = 1:length(detectMonth)
    for ii = 1:size(detectMonth,1)
        if ~isempty(detectMonth{ii,i})
            for iii = 1:length(detectMonth{ii,i})
                if ~isempty(detectMonth{ii,i}(iii).levelVec)
                    date = [num2str(month(detectMonth{ii,i}(iii).timeStamp(1))), '_', num2str(year(detectMonth{ii,i}(iii).timeStamp(1)))];
                    name = strcat(date,'_detection');
                    fileName = fullfile(dataLoc,boxNames{i},compDataLoc,name);
                    levelVec = detectMonth{ii,i}(iii).levelVec;
                    totalData = detectMonth{ii,i}(iii).totalData;
                    responseVec = detectMonth{ii,i}(iii).responseVec;
                    timeStamp = detectMonth{ii,i}(iii).timeStamp;
                    toneVec = detectMonth{ii,i}(iii).toneVec;
                    phaseChoice = detectMonth{ii,i}(iii).phaseChoice;
                    target = detectMonth{ii,i}(iii).target;
                    nontarget = detectMonth{ii,i}(iii).nontarget;
                    save(fileName,'levelVec','totalData','responseVec','timeStamp','toneVec','phaseChoice','target','nontarget')
                    clear date name fileName levelVec totalData responseVec timeStamp toneVec phaseChoice target nontarget
                else
                end
            end
        end
    end
end
%save monthly discrimination data into separate files by tonebox
for i = 1:length(discrimMonth)
    if ~isempty(discrimMonth{i})
        for ii = 1:length(discrimMonth{i})
            date = [num2str(month(discrimMonth{i}(ii).timeStamp(1))), '_', num2str(year(discrimMonth{i}(ii).timeStamp(1)))];
            name = strcat(date,'_discrimination');
            fileName = fullfile(dataLoc,boxNames{i},compDataLoc,name);
            levelVec = discrimMonth{i}(ii).levelVec;
            totalData = discrimMonth{i}(ii).totalData;
            responseVec = discrimMonth{i}(ii).responseVec;
            timeStamp = discrimMonth{i}(ii).timeStamp;
            toneVec = discrimMonth{i}(ii).toneVec;
            phaseChoice = discrimMonth{i}(ii).phaseChoice;
            target = discrimMonth{i}(ii).target;
            nontarget = discrimMonth{i}(ii).nontarget;
            save(fileName,'levelVec','totalData','responseVec','timeStamp','toneVec','phaseChoice','target','nontarget')
            clear date name fileName levelVec totalData responseVec timeStamp toneVec phaseChoice target nontarget
        end
    end
end
