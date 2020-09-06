%% Script to post-process dF/F traces from volume data
% Functions at the bottom may make some assumptions that only apply to tone
% experiments produced via Psignal.
% Zac Bowen

clear
datadir =  {'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031020_367n_100um20st_FRA\' % 1
            'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031020_367r_100um20st_FRA\' % 2
            'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031020_ia09_100um20st_FRAinj\' % 3
            'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031120_352ll_100um20st_FRA\' % 4
            'C:\Users\Zac\Dropbox (Personal)\research\volume_data\031120_352ll_100um20st_FRA_diffxy\'}; % 5
        
%% create processed data files
for fileNum = 1:length(datadir)
    clearvars -except datadir fileNum
    cd(datadir{fileNum})

    XML = danParseXML('Experiment.xml');
    exptVars = xmlVersionCheck(XML,[]);
    [zDFF,allxc,allyc,allzc] = deal(cell(exptVars.totalZplanes,1));
    for z = 1:exptVars.totalZplanes
        load(['DFF_Z' num2str(z) '.mat'],'DFF','xc','yc','keptRois')
        % Typically 
        zDFF{z} = DFF.raw;
        clear DFF %RASTERS
        
        % I've stopped using the "cleaned" traces since I updated my computeDFF function
%         xc = xcRaw(keptRois.rawClean); 
%         yc = ycRaw(keptRois.rawClean);
        
        % Compute z-coords
        stepSizePix = exptVars.stepSizeUM/exptVars.micronsPerPixel;
        zFun = (stepSizePix:-(stepSizePix/512):(stepSizePix/512))-(z*stepSizePix);
        zc = zFun(round(yc))';
        
        allxc{z} = xc;
        allyc{z} = yc;
        allzc{z} = zc;
        allkeptRois(z) = keptRois;
        [BFinfo(z),CorrInfo(z),CellInfo(z),zStuff(z)] = getTuningPsignalZ(z,exptVars.totalZplanes);
    end

    [allZCorrInfo] = getCorrsAllZ(zStuff);
    selectedZplanes = 2:6; % *** Other users may want to change this
    [selectZCorrInfo] = getCorrsAllZ(zStuff(selectedZplanes));

    save(['allPlanesVariables' date '.mat'],...
        'allZCorrInfo','selectZCorrInfo','zDFF','allxc','allyc','allzc',...
        'BFinfo','CorrInfo','CellInfo','zStuff','exptVars')
end

%% Functions
function [BFinfo,CorrInfo,CellInfo,zStuff] = getTuningPsignalZ(zNum,numZplanes)


filename = ['DFF_Z' num2str(zNum) '.mat'];
if ~isfile(filename)
    fprintf('File not in directory\n')
    BFinfo = [];
    return
end
fprintf('loading %s...',filename)
load(filename)
fprintf('done\n')

psigfile = dir('RND*.mat');
stimInfo =  WF_getPsignalInfo(fullfile(psigfile.folder,psigfile.name));
load('trialStartFrames.mat')
trialStartVals = round(trialStartFrame/numZplanes);
trialEndVals = trialStartVals + round(median(diff(trialStartFrame))/numZplanes);

ImgFrameRate = 30; % in Hz (frames/sec) -- Assume for these datasets
planeFrameRate = ImgFrameRate/numZplanes;
analysisWinSeconds = 1; % in Seconds
analysisWinFrames = round(analysisWinSeconds*planeFrameRate);


nstims = [numel(stimInfo.uFreqs) numel(stimInfo.uLevels)];
kHzVals = stimInfo.uFreqs/1000;
repOrdIdx = stimInfo.Trialindicies(:,1);
% stimFrame = round(stimInfo.PreStimSilence*planeFrameRate);
stimFrame = round(trialStartVals + stimInfo.PreStimSilence*planeFrameRate);
% stimFrame = round(trialStartFrame/numZplanes + stimInfo.PreStimSilence*planeFrameRate); % possibly do this
% dff = DFF.smoothClean; % dff should be ROIs x Frames
% xc = xc; yc = yc;
% xc = xcRaw; yc = ycRaw;

% dff = DFF.rawClean; % dff should be ROIs x Frames
% xc = xcRaw(keptRois.rawClean); yc = ycRaw(keptRois.rawClean);
dff = DFF.raw; % dff should be ROIs x Frames

numStims = length(unique(repOrdIdx));
numTrials = length(repOrdIdx)/numStims;
numCells = size(dff,1);

pTrial = NaN(numCells,numStims,numTrials);
pStim = NaN(numCells,numStims);
pOff = NaN(numCells,numStims);
for stim = 1:numStims
    trials = find(repOrdIdx==stim);
    for tr = 1:length(trials)
        trialDFF(stim,tr,:,:) = dff(:,trialStartVals(trials(tr)):trialEndVals(trials(tr)));
        baseFrames = (stimFrame(trials(tr)) - 1 - analysisWinFrames) : (stimFrame(trials(tr)) - 1); %30:60 if stimFrame=61
        stimFrames = (stimFrame(trials(tr)))                         : (stimFrame(trials(tr)) + analysisWinFrames); %61:91
        offFrames  = (stimFrame(trials(tr)) + analysisWinFrames + 1) : (stimFrame(trials(tr)) + 2*analysisWinFrames + 1); %61:91
%         if offFrames(end) > range(trialStartVals(trials(tr)):trialEndVals(trials(tr))); fprintf('Warning: offFrames go into next trial\n'); end
        baseDFF(stim,tr,:,:) = dff(:,baseFrames);
        stimDFF(stim,tr,:,:) = dff(:,stimFrames);
        offDFF(stim,tr,:,:)  = dff(:,offFrames);
%         baseDFF(stim,tr,:,:) = trialDFF(stim,tr,:,baseFrames);
%         stimDFF(stim,tr,:,:) = trialDFF(stim,tr,:,stimFrames);
%         offDFF(stim,tr,:,:)  = trialDFF(stim,tr,:,offFrames);
        
        for c = 1:numCells
            if (mean(baseDFF(stim,tr,c,:)) < mean(stimDFF(stim,tr,c,:)))% && (mean(stimDFF(stim,tr,c,:)) > 3) % at least 3% dff for response
                [pTrial(c,stim,tr),~] = statTest2(squeeze(baseDFF(stim,tr,c,:)), squeeze(stimDFF(stim,tr,c,:)));
            end
        end
        
    end
    
    for c = 1:numCells
        trialMean_baseDFF = mean(baseDFF(stim,:,c,:),4)'; % Tried concatenating instead of mean
        trialMean_stimDFF = mean(stimDFF(stim,:,c,:),4)'; % way too many stims were significant (even at 0.05)
        trialMean_offDFF = mean(offDFF(stim,:,c,:),4)';
        if (mean(trialMean_baseDFF) < mean(trialMean_stimDFF))
            [pStim(c,stim),~] = statTest2(trialMean_baseDFF, trialMean_stimDFF);
        end
        if ~(pStim(c,stim)<0.05) && (mean(trialMean_baseDFF) < mean(trialMean_offDFF))
            [pOff(c,stim),~] = statTest2(trialMean_baseDFF, trialMean_offDFF);
        end
        flatFRA(c,stim) = mean(trialMean_stimDFF); % cells x unique stimuli (vector)
        flatTrialFRA(c,stim,:) = trialMean_stimDFF; % cells x unique stim x trial mean
    end
    
end
sigTrial = zeros(numCells,numStims,numTrials);
sigTrial(pTrial<0.05) = 1;
sigStim = zeros(numCells,numStims);
sigStim(pStim<0.05) = 1; % 0.05 or 0.01?
sigOff = zeros(numCells,numStims);
sigOff(pOff<0.05) = 1; % 0.05 or 0.01?

a = reshape(flatFRA,size(flatFRA,1),nstims(1),nstims(2));
FRA = permute(a,[1 3 2]);
aa = reshape(sigStim,size(sigStim,1),nstims(1),nstims(2));
sigStimFRA = permute(aa,[1 3 2]);

%compare this sigStim result with if you used frame/trial concatenation
%rather than mean of frames. or increase the pvalue for freqs
%%
[BFval,BL,BFresp,CFresp,CFval,bandwidth] = deal(NaN(numCells,1));
[selectivity,selectivityBinary] = deal(NaN(numCells,size(FRA,2)));
for c = 1:numCells

    unsigInds = find(~squeeze(sigStimFRA(c,:,:))); %freqs that had no sig response
    
    % Calculate selectivity measure
    tmp = squeeze(FRA(c,:,:));
    tmp(unsigInds) = 0; %set non-signif resps to 0
    %had to make this weird if-statement for all-negative FRAs.
    %prob a smarter way to do all of this but it works
    if max(tmp(:))<=0; tmp=abs(tmp+min(tmp(:))); end
    tmp2 = tmp - min(tmp(:));
    tmpnormFRA = tmp2/max(tmp2(:)); %normalize FRA of just sig responses
    tmpnormFRAbinary = tmpnormFRA > 0; % binarized FRA
    for L = 1:size(tmpnormFRA,1)
        selectivity(c,L) = sum(tmpnormFRA(L,:)); % sum the normalized significant responses
        selectivityBinary(c,L) = sum(tmpnormFRAbinary(L,:)); % sum the normalized significant responses
%     selectivity(c) = sum(tmpnormFRA(:)); % sum the normalized significant responses
    end
    
    % Calculate characteristic frequency (CF)
    lowestSigLevel = find(sum(tmp,2),1,'last'); % finds the lowest dB Level with sig resp
    if ~isempty(lowestSigLevel)
        [CFresp(c) , CFval(c) ] = max( tmp(lowestSigLevel,:) );
    end
    
    % Calculte best frequency (BF)
    [tmpResp , tmpBF ] = max( flatFRA(c,:) );
    BFval(c) = mod(tmpBF,nstims(1));  % Preallocate
    if mod(tmpBF,nstims(1))==0
        BFval(c) = nstims(1);
    end
    BL(c) = ceil(tmpBF/nstims(1));
    BFresp(c) = tmpResp;
    
    % Save normalized FRA (probably not necessary)
    tmp = squeeze(FRA(c,:,:));
    tmp2 = tmp - min(tmp(:));
    tmpnormFRA = tmp2/max(tmp2(:));
    normFRA(c,:,:) = tmpnormFRA;

% the first way I tried to do the selectivity. modified to Ying's way above
%     sigInds = find(squeeze(sigStimFRA(c,:,:)));
%     onlySigResps = tmpnormFRA(sigInds);
%     selectivity(c) = sum(onlySigResps);
    if any(isnan(squeeze(FRA(c,BL(c),:))))
        bandwidth(c) = NaN;
    else
        bandwidth(c) = peak_analyzer(squeeze(FRA(c,BL(c),:)),0.7,kHzVals);
    end
end

respCells = find(sum(sigStim,2)~=0);
unrespCells = find(sum(sigStim,2)==0);

%% Calculate signal and noise correlations
nTrials = size(flatTrialFRA,3);
SigCorrs = NaN(numCells,numCells);
NoiseCorrsTrial = NaN(numCells,numCells,nTrials);
NoiseCorrsVec = NaN(numCells,numCells);
pairs = nchoosek(1:numCells,2); % Define all unique pairwise relationships
for i = 1:size(pairs,1) % Loop through pairs
    
    c1 = pairs(i,1); % Define each cell
    c2 = pairs(i,2);
    TC_1 = flatFRA(c1,:); % Define FRA of each cell
    TC_2 = flatFRA(c2,:);
    C = cov(TC_1,TC_2); % Covariance of the two FRAs
    P_coeff = C(2)/sqrt(C(1)*C(4)); % Correlation coefficient
    % P_coeff = C(2)/(std(TC_1)*std(TC_2)); % Alternate form of calculation 
    SigCorrs(c1,c2) = P_coeff;
    
    for tr = 1:nTrials % Loop through trials
        % Do same as above but use trial FRA deviations for the cross correlation
        TC_1_nc = flatFRA(c1,:) - flatTrialFRA(c1,:,tr);
        TC_2_nc = flatFRA(c2,:) - flatTrialFRA(c2,:,tr);
        C_nc = cov(TC_1_nc,TC_2_nc);
        P_coeff_nc = C_nc(2)/sqrt(C_nc(1)*C_nc(4));
        NoiseCorrsTrial(c1,c2,tr) = P_coeff_nc;
    end
    
    TC_1_ncAlt = reshape(flatTrialFRA(c1,:,:),1,nTrials*numStims) ...
                - repmat(flatFRA(c1,:),1,nTrials);
    TC_2_ncAlt = reshape(flatTrialFRA(c2,:,:),1,nTrials*numStims) ...
                - repmat(flatFRA(c2,:),1,nTrials);
    C_ncAlt = cov(TC_1_ncAlt,TC_2_ncAlt);
    P_coeff_ncAlt = C_ncAlt(2)/sqrt(C_ncAlt(1)*C_ncAlt(4));
    NoiseCorrsVec(c1,c2) = P_coeff_ncAlt;
end

%% Save relative pairwise cell locations
cellDists = pairwiseDist(xc,yc); % in pixels
cellAngles = pairwiseAngle(xc,yc); % in degrees

%% Organize all of the results
BFinfo.cellnum = respCells;
BFinfo.NONcellnum = unrespCells;
BFinfo.BFval = BFval;
BFinfo.BFresp = BFresp;
BFinfo.CFval = CFval;
BFinfo.CFresp = CFresp;
BFinfo.BL = BL;
BFinfo.fraVals = FRA;
BFinfo.bandwidth = bandwidth;
BFinfo.normFRA = normFRA;
BFinfo.selectivity = selectivity;
BFinfo.selectivityBinary = selectivityBinary;
BFinfo.sigRespCells = find((sum(sigStim,2)));
BFinfo.sigOffRespCells = find((sum(sigOff,2)));
CorrInfo.SigCorrs = SigCorrs;
CorrInfo.NoiseCorrsTrial = NoiseCorrsTrial;
CorrInfo.NoiseCorrsVec = NoiseCorrsVec;
CellInfo.cellDists = cellDists;
CellInfo.cellAngles = cellAngles;
CellInfo.sigTrial = sigTrial;
CellInfo.sigStim = sigStim;
CellInfo.sigOff = sigOff;
zStuff.flatFRA = flatFRA;
zStuff.flatTrialFRA = flatTrialFRA;
zStuff.pStim = pStim;
zStuff.pTrial = pTrial;

save(['tuningVarsZ' num2str(zNum) '.mat'],'nstims','trialDFF','stimFrame','trialStartVals','FRA','kHzVals')
% save(filename,'BFinfo','CorrInfo','CellInfo','-append')
end


function [CorrInfo] = getCorrsAllZ(zStuff)

flatFRA = [];
flatTrialFRA = [];
for z = 1:length(zStuff)
    flatFRA      = cat(1,flatFRA     ,zStuff(z).flatFRA);
    flatTrialFRA = cat(1,flatTrialFRA,zStuff(z).flatTrialFRA);
end
numCells = size(flatFRA,1);

%% Calculate signal and noise correlations
nTrials = size(flatTrialFRA,3);
SigCorrs = NaN(numCells,numCells);
NoiseCorrsTrial = NaN(numCells,numCells,nTrials);
% NoiseCorrsVec = NaN(numCells,numCells);
pairs = nchoosek(1:numCells,2); % Define all unique pairwise relationships
for i = 1:size(pairs,1) % Loop through pairs
    
    c1 = pairs(i,1); % Define each cell
    c2 = pairs(i,2);
    TC_1 = flatFRA(c1,:); % Define FRA of each cell
    TC_2 = flatFRA(c2,:);
    C = cov(TC_1,TC_2); % Covariance of the two FRAs
    P_coeff = C(2)/sqrt(C(1)*C(4)); % Correlation coefficient
    % P_coeff = C(2)/(std(TC_1)*std(TC_2)); % Alternate form of calculation 
    SigCorrs(c1,c2) = P_coeff;
    
    for tr = 1:nTrials % Loop through trials
        % Do same as above but use trial FRA deviations for the cross correlation
        TC_1_nc = flatFRA(c1,:) - flatTrialFRA(c1,:,tr);
        TC_2_nc = flatFRA(c2,:) - flatTrialFRA(c2,:,tr);
        C_nc = cov(TC_1_nc,TC_2_nc);
        P_coeff_nc = C_nc(2)/sqrt(C_nc(1)*C_nc(4));
        NoiseCorrsTrial(c1,c2,tr) = P_coeff_nc;
    end
    
%     TC_1_ncAlt = reshape(flatTrialFRA(c1,:,:),1,nTrials*numStims) ...
%                 - repmat(flatFRA(c1,:),1,nTrials);
%     TC_2_ncAlt = reshape(flatTrialFRA(c2,:,:),1,nTrials*numStims) ...
%                 - repmat(flatFRA(c2,:),1,nTrials);
%     C_ncAlt = cov(TC_1_ncAlt,TC_2_ncAlt);
%     P_coeff_ncAlt = C_ncAlt(2)/sqrt(C_ncAlt(1)*C_ncAlt(4));
%     NoiseCorrsVec(c1,c2) = P_coeff_ncAlt;
end

%% Organize all of the results
CorrInfo.SigCorrs = SigCorrs;
CorrInfo.NoiseCorrsTrial = NoiseCorrsTrial;
% CorrInfo.NoiseCorrsVec = NoiseCorrsVec;

% save(filename,'BFinfo','CorrInfo','CellInfo','-append')
end