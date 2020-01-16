%This script creates analyzes tonebox detection-phase behavior that has
%already been split up according to individual dates for each tonebox.
%MUST BE RUN AFTER "toneBoxBehaviorDailySetup" TO ENSURE DAILY DATA FILES
%EXIST IN THE CORRECT LOCATION!
%Ilan Goldstein (11/2019)
box = input('Select tonebox for analysis: ');
tonebox = strcat('tonebox',num2str(box));
genPath = 'C:\Users\PsiDev\Desktop\Ilan_Psignal';                          %location of relevant psignal code and this script
dataLoc = 'D:\Devices';                                                    %parent folder containing tonebox data folders
DDloc = fullfile(dataLoc,tonebox,'daily_data');                            %location of daily tonebox data (stored by "toneBoxBehaviorDailySetup")
cd(DDloc)
DDfiles = dir('*detection*');                                              %list of all daily behavior data files
%%load each daily detection data file and run through%%
%%"toneBoxAnalysisMonthly" to output relevant behavior data%%
for i = 1:length(DDfiles)
    DDload = DDfiles(i).name;
    DDloadFile = fullfile(DDloc,DDload);
    [bhvOutput] = toneBoxAnalysisMonthly(DDloadFile);                      %"toneBoxAnalysisMonthly" (TBAM) can be used to analyze monthly or daily data files
    set(gcf, 'WindowStyle', 'Docked')
    date(i) = datetime(extractBefore(DDload,12));                          %combining TBAM output data into single-tonebox matrices
    earlyRates(i) = bhvOutput.earlyRate;
    hitRates(i) = bhvOutput.hitRate;
    falseAlarmRates(i) = bhvOutput.falseAlarmRate;
    dprimes(i) = bhvOutput.dprime;
    targetResps(i,:) = bhvOutput.targetResp(:,2);
    nontargetResps(i,:) = bhvOutput.nontargetResp(:,2);
    t = bhvOutput.t;
    lickResponseE(i,:) = bhvOutput.lickResponseE;
    lickResponseH(i,:) = bhvOutput.lickResponseH;
    lickResponseF(i,:) = bhvOutput.lickResponseF;
end
%%sort data matrices according to sorted dates of daily data files%%
[srtd didx] = sort(date);
sDate = date(didx);
sER = earlyRates(didx)';
sHR = hitRates(didx)';
sFR = falseAlarmRates(didx)';
sDprimes = dprimes(didx);
sTarResps = targetResps(didx,:);
sNonResps = nontargetResps(didx,:);
sLRE = lickResponseE(didx,:);
sLRH = lickResponseH(didx,:);
sLRF = lickResponseF(didx,:);
respRates = [sER sHR sFR];
%%calculate date value (start date = 1) for all daily data files%%
for i = 1:length(sDate)
    dur = between(sDate(1),sDate(i),'Days');
    diff(i) = str2num(extractBefore(char(dur),'d')) + 1;
end
%%plot average daily data across all detection dates for the tonebox being analyzed%%
figure
suptitle([tonebox,' daily detection'])                                     %t: time intervals used by data collection software to match trial length              
subplot(1,8,1)                                                             %diff: date value vector for all dates being analyzed
title('avg target response')
waterfall(t,diff,sTarResps)                                                %daily average target-tone-trial lick responses
subplot(1,8,2)
title('avg nontarget response')
waterfall(t,diff,sNonResps)                                                %daily average nontarget-tone-trial lick responses
subplot(1,8,3)
title('avg response latency: early')
waterfall(t,diff,sLRE)                                                     %daily average early response latencies
subplot(1,8,4)
title('avg response latency: hit')
waterfall(t,diff,sLRH)                                                     %daily average hit response latencies
subplot(1,8,5)
title('avg response latency: false alarm')
waterfall(t,diff,sLRF)                                                     %daily average false alarm response latencies
subplot(1,8,6:8)
bar(respRates)                                                             %daily average response rates
set(gcf, 'WindowStyle', 'Docked')