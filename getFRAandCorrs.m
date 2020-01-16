
% This function loads the .mat datafile (after dF/F calculation) specified by 
% datadir/exp/filename, computes FRAs and corrs, then appends the .mat file
function BFinfo = getFRAandCorrs(input)

%% --- Front end crap begins ---
 
% fps = input.expectedFPS;
% for e = 1:length(input.expname)
%     PsignalFile = input.psignalfiles{e};
%     bb = strsplit(input.path,'\');
%     InPath   = fullfile(input.path , input.expname{e}) ;
%     SavePath = input.savepath;
%     SavePath = [SavePath '\'  bb{5}  '\' input.expname{e} '\'];
%    % copyfile ( InPath, SavePath);
%     PsignalFile=[SavePath PsignalFile];
%     handles = WF_getPsignalInfo(PsignalFile);
%     try
%     Output = load(fullfile( InPath, 'Fluorescence.mat'));
%     catch 
%         warning( 'No Fluorescence Data!')
%         continue
%     end 
%     FCellCorrected = Output.FCellCorrected;
%     trialdur = size(FCellCorrected,1);
%     Neurons  = size(FCellCorrected,3);
% 
%  
%      % TO-DO Fix the Creation of the psignal matrix
%     L      = length(FCellCorrected);
%     uFreqs  = handles.uFreqs ;
%     uF    = length(uFreqs);
%     uLevels = handles.uLevels;
%     freqs   = handles.Freqs  ;
%     freqs   = freqs(1:L);
%     Levels  = handles.Levels ;
%     Levels  = Levels(1:L);
%     FreqLevelOrder = table(freqs,Levels);
%     FreqLevels = unique(FreqLevelOrder);
% 
% 
% 
% 
% fprintf(['Loading ' filename '...'])
        
% %         stimFrame =  % Frame of stim onset
%         stimStart = stimFrame / ImgFrameRate; % Stim onset in seconds
%         gateDur = handles.PrestimSilence + handles.PostStimSilence+...
%                   handles.PrimaryDuration; % Trial duration in seconds
%         freqs = uF/1000),1); % Stim values in kHz
%         levels = unique(TagMatrix(stimFrame,:,strcmp(Tags,'ToneLevel'))); % Stim levels in dB
%         nstims = [numel(freqs) numel(levels)];
%         freqorder = TagMatrix(stimFrame,:,strcmp(Tags,'ToneFrequency'));
%         repOrdIdx = round(2*(log2(freqorder/min(freqorder(:)))+0.5));
% end

%% --- Front end crap ends ---

% inputDFF should be nFrames x nROIs
% inputDFF = DFF.smoothClean';
% nROIs = size(inputDFF,2);
% 
% % "nstims" is 2 element array of [numLevels numFreqs]
% numStims = nstims(1)*nstims(2);
% 
% % parsedFluoTraces is ( total stims x trial length x nROIs)
% parsedDFF = cell(1, numStims);% , tarFnums+1 , nROIs );

%% Parse the individual neuron traces into trial based events/cells
for stimNum = 1:(nstims(1)*nstims(2)) % Loop through unique stimuli
    
    traceIdx = find(repOrdIdx == stimNum); % Find presentations of this stimuli
    
    for rCtr = 1:length(traceIdx) % Loop through those trials and pull out dF/F traces

        parsedDFF{stimNum} (rCtr , : , : ) = inputDFF ( seqStartVals ( traceIdx(rCtr) ) : seqStartVals ( traceIdx(rCtr) )+tarFnums , : );

    end
end


%% SIGNFICANCE TESTING - DW08122015
%% DOES THE NEURON RESPOND?

% Define Analysis Window length for stim frames and base frames
AnalysisWin = 1; % in Seconds
analysisWinFrames = AnalysisWin * ImgFrameRate;

sigFlag = zeros(1,nROIs);
tunFlag = zeros(1,nROIs);
% sigFlag indicates if a neuron significantly responds to ANY stimulus.
% tunFlag indicates if a neuron was tuned (specific in its response)

for ii= 1:nROIs % Loops through all ROIs
    tmpROI=[];
    for kk = 1:(nstims(1)*nstims(2))
        tmpROI = [tmpROI;parsedDFF{kk}(:,:,ii)];
    end
    %%% ZCB 10/27/16 started to debug 1 sec prestim here
    %%% Problem is that you're using 29 frames instead of the usual 30
    %%% Your prestim period ends up starting at index 0. matlab starts at 1
    preFrames_start = median(stimFrame) - analysisWinFrames;
    if preFrames_start == 0
        preFrames_start = 1;
    end
    %     preStim = tmpROI(:,preFrames_start : median(stimFrame)-1);%grabPreStimulusFrames % 11/16/16 Commented out for line below for new data structure
    preStim = tmpROI(:,median(stimFrame)-analysisWinFrames + 1  : median(stimFrame));%grabPreStimulusFrames
    postStim = tmpROI(:,median(stimFrame)+1 : median(stimFrame)+analysisWinFrames);%grabPostStimulusFrames
    
    if mean(postStim(:)) > 3 %hard-coded threshold for defining a response - at least 3% dff
        [p,~,stats] = anova1( [preStim(:) postStim(:)] , [] , 'off' );
        if p < 0.01 % DOES IT PREFER ANY STIMULUS - NOT WHICH STIMULUS
            sigFlag(ii) = 1;
            
            clear tmpTun
            for zz = 1:(nstims(1)*nstims(2))
                tmpTun(zz,:,:,:) = parsedDFF{zz}(:,:,ii);
            end
            
            [pTun,~,stats] = ...
                anova1( ... % , [] , 'off' );
                  mean( ... % ,3)'
               squeeze( ... % )
                        tmpTun(:,:,median(stimFrame)+1:median(stimFrame)+analysisWinFrames) ...
                        ),3)',[],'off');
                    
            if pTun < 0.01
                tunFlag(ii) = 1;
            end
            
        end
    end
    
end


%% FIND THE BFs FOR EACH NEURON THAT RESPONDS AND IS TUNED - MID-RANGE INTENSITIES
   idx = find(  sigFlag == 1 & tunFlag==1 ); % find the neurons with sigResp and sigTuning
NONidx = find(~(sigFlag == 1 & tunFlag==1));
for nn = 1:length(idx) % Loops through tuned and selective cells
    clear tmpTun
    for zz = 1:(nstims(1)*nstims(2))
        tmpTun(zz,:,:,:) = parsedDFF{zz}(:,:,idx(nn));
    end
    
    fraRange = (nstims(1)*(nstims(2)-1))-(nstims(1)-1) : (nstims(1)*(nstims(2)-1));
    
    [tmpResp , tmpBF ] = ...
        max( ... % )
        mean( ... %, 2)
        mean( ... %, 3)
             tmpTun( fraRange, :, (median(stimFrame)+1 : median(stimFrame)+analysisWinFrames) ) ...
             ,3),2));
    
    BFval(nn) = tmpBF;
    BFresp(nn) = tmpResp;
        
end 

%% Calculate average responses for each stimulus to be used in FRA
for yy = 1:numStims
    % Average response across stim frames averaged over trials
    avgResp(yy,:) = ...
        squeeze( ... % )
           mean( ... %, 2)  across frames
           mean( ... %, 1)  across reps
                parsedDFF{yy} ( : , (median(stimFrame)+1) : (median(stimFrame)+ (AnalysisWin * ImgFrameRate)), :) ...
                ,2) ,1) );
            
    % Average response across stim frames seperated for each trial
    avgTrialResp(yy,:,:) = ...
        squeeze( ... % )
           mean( ... %, 2)  across frames
                parsedDFF{yy} ( : , (median(stimFrame)+1) : (median(stimFrame)+ (AnalysisWin * ImgFrameRate)), :) ...
                ,2) );
            
%     avgTCResp(yy,:,:) = ...
%         squeeze( ... % )
%          median( ... %, 1)
%                 parsedDFF{yy}(:,:,:) ...
%                 ,1));

end

%% Organize FRAs from the averaged responses (cell ID, db Lvl, freq)
for xx = 1:nROIs % Loops through cells
    fraVals(xx,:,:) = reshape(avgResp(:,xx), nstims(1), nstims(2))';
    for tr = 1:size(avgTrialResp,2)
        fraTrialVals(xx,:,:,tr) = reshape(avgTrialResp(:,tr,xx), nstims(1), nstims(2))';
    end
    %     bandwidth(xx) = peak_analyzer(squeeze(fraVals(xx,2,:)),0.7,freqs); % custom function by ZCB
    normFRA(xx,:,:) = fraVals(xx,:,:) - min(min(fraVals(xx,:,:)));
    normFRA(xx,:,:) = normFRA(xx,:,:) / max(max(normFRA(xx,:,:)));
end

%% Assign desired tuning info to output structure "BFinfo"
BFinfo.cellnum = idx;
BFinfo.NONcellnum = NONidx;
BFinfo.BFval = BFval;
BFinfo.BFresp = BFresp;
BFinfo.fraVals = fraVals;
BFinfo.bandwidth = bandwidth;
BFinfo.normFRA = normFRA;

%% Compute signal and noise correlations using tuning curves -- ZCB added late 2017
originalFRAs = fraVals;
tmpfraT = permute(originalFRAs,[1 3 2]);
originalFRAsflat = reshape(tmpfraT,size(tmpfraT,1),size(tmpfraT,2)*size(tmpfraT,3));

originalTrialFRAs = fraTrialVals;
tmpTrialfraT = permute(originalTrialFRAs,[1 3 2 4]);
originalTrialFRAsflat = reshape(tmpTrialfraT,size(tmpTrialfraT,1),size(tmpTrialfraT,2)*size(tmpTrialfraT,3),size(tmpTrialfraT,4));

nROIs = size(originalFRAs,1);
nTrials = size(originalTrialFRAsflat,3);
SigCorrs = NaN(nROIs,nROIs);
NoiseCorrsTrial = NaN(nROIs,nROIs,nTrials);
NoiseCorrsVec = NaN(nROIs,nROIs);
pairs = nchoosek(1:nROIs,2); % Define all unique pairwise relationships
for i = 1:size(pairs,1) % Loop through pairs
    
    c1 = pairs(i,1); % Define each cell in pair
    c2 = pairs(i,2);
    
    % Signal correlations
    TC_1 = originalFRAsflat(c1,:); % Define FRA of each cell
    TC_2 = originalFRAsflat(c2,:);
    C = cov(TC_1,TC_2); % Covariance of the two FRAs
    P_coeff = C(2)/sqrt(C(1)*C(4)); % Correlation coefficient
    % P_coeff = C(2)/(std(TC_1)*std(TC_2)); % Alternate form of calculation 
    SigCorrs(c1,c2) = P_coeff;
    
    % Noise correlations on trial-by-trial FRA deviations within each stimulus "block"
    % This gives you a noise corr value for each stim block (aka number of
    % trial repeats)
    for tr = 1:nTrials % Loop through trials
        % Do same as above but use trial FRA deviations for the cross correlation
        TC_1_nc = originalFRAsflat(c1,:) - originalTrialFRAsflat(c1,:,tr);
        TC_2_nc = originalFRAsflat(c2,:) - originalTrialFRAsflat(c2,:,tr);
        C_nc = cov(TC_1_nc,TC_2_nc);
        P_coeff_nc = C_nc(2)/sqrt(C_nc(1)*C_nc(4));
        NoiseCorrsTrial(c1,c2,tr) = P_coeff_nc;
    end
    
    % Noise correlations on trial-by-trial FRA deviations vectorized
    % This gives you one noise corr value since all trial blocks are
    % concatenated prior to calculation (ZCB preference)
    TC_1_ncAlt = reshape(originalTrialFRAsflat(c1,:,:),1,nTrials*numStims) ...
                - repmat(originalFRAsflat(c1,:),1,nTrials);
    TC_2_ncAlt = reshape(originalTrialFRAsflat(c2,:,:),1,nTrials*numStims) ...
                - repmat(originalFRAsflat(c2,:),1,nTrials);
    C_ncAlt = cov(TC_1_ncAlt,TC_2_ncAlt);
    P_coeff_ncAlt = C_ncAlt(2)/sqrt(C_ncAlt(1)*C_ncAlt(4));
    NoiseCorrsVec(c1,c2) = P_coeff_ncAlt;
    
    % Noise correlations could also be computed from a covariance of
    % spontaneous portions of dF/F. Ideally you would want a significant
    % spontaneous "block" either at the end or beginning of your experiment.
    % i.e. Dan would have ~8min spontaneous at the end of stim experiments
    
end

% How do trial-block-specific NCs compare to a mean across trials?
% figure;
% subplot(ceil(sqrt(nTrials)),floor(sqrt(nTrials)),1)
% histogram(nanmean(NoiseCorrsTrial,3))
% axis('square'); title('Mean across trials')
% for tr = 1:nTrials
%     subplot(ceil(sqrt(nTrials)),floor(sqrt(nTrials)),tr+1)
%     histogram(NoiseCorrsTrial(:,:,tr))
%     axis('square'); title(['trial ' num2str(tr)])
% end

f.CorrInfo.SigCorrs = SigCorrs;
f.CorrInfo.NoiseCorrsTrial = NoiseCorrsTrial;
f.CorrInfo.NoiseCorrsVec = NoiseCorrsVec;
f.BFinfo = BFinfo;

% Uncomment these lines for saving -- careful not to overwrite your
% datafiles while you're testing out this code
% % % % % Save the old file as backup before appending new variables
% % % % if isempty(dir(['backup' date filename]))
% % % %     copyfile(filename,['backup' date filename])
% % % % end
% % % % 
% % % % % save(filename,'BFinfo','CorrInfo','-append')
% % % % save(filename,'f','-append')


end