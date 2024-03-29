function ExtractDeltaFSpikes(input)

for e = 1:length(input.expname)
    PsignalFile = input.psignalfiles{e};
    bb = strsplit(input.path,'\');
    ExpName = input.expname{e};
    LoadPath = [input.savepath '\'  bb{5}  '\' ExpName '\'];
    PsignalPath = [LoadPath];
    disp(['Extracting DeltaF and Deconvolving spikes for ' PsignalFile])
    ThorFile = [LoadPath '\timing.txt'];
    %Load data
    load([LoadPath 'PsignalMatrix.mat'])
    load([LoadPath 'Fluorescence.mat'])
    %Extract timing parameters
    TimingInfo = getTimingInfo(ThorFile, [PsignalPath PsignalFile], input.expectedFPS);
    %Find each trial's tone onset frame and frequency
    TagNames = categorical(PsignalMatrix.TagNames(:,2));
    idx = find(TagNames=='StimOnset');
    Freqs = squeeze(PsignalMatrix.Tags(:,:,idx));
    F = unique(Freqs(:));
    F = F(F>0);
    toneonsetframe = [];
    ftrial=[];
    for i = 1:size(Freqs,2)
        toneonsetframe(i) = find(Freqs(:,i));
        ftrial(i) = Freqs(toneonsetframe(i),i);
    end
    %Find cell brightness (i.e., cell ring flouro re. neuropil)
    FCell = Output.FCell;
    FNeuropil = Output.FNeuropil;
    cellbrightness=squeeze(100*(nanmean(nanmean(FCell,2),1)-nanmean(nanmean(FNeuropil,2),1))./nanmean(nanmean(FNeuropil,2),1));
    brightcells = cellbrightness>0;
    %Smooth fluorescence using Plenz lab approach
    T = 1/input.expectedFPS;
    win = 10; %Change based on fast rise and exponential decay
    tau = 1.5;
    wf = [0:1:win];
    wf = exp(-tau*T*wf);
    FCellCorrectedSmooth=[];
    FCellCorrected = Output.FCellCorrected;
    for i = 1:size(FCellCorrected,3)
        for ii = 1:size(FCellCorrected,2)
            tmp = FCellCorrected(:,ii,i)';
            % 1 ->symmetric, 2->asym, T = frame time, Tau -> 1.5
            symmFLAG = 1;
            % Because T == 0.033 always
            smoothwin = 9;
            FCellCorrectedSmooth(:,ii,i) = smooth2(tmp,symmFLAG,smoothwin,T,1.5);
        end
    end
    Output.FCellCorrectedSmooth = FCellCorrectedSmooth;
    %Find DeltaF
    DeltaFCellCorrected=[];
    DeltaFCellCorrectedSmooth=[];
    for i = 1:size(Freqs,2)
        baseline = nanmean(squeeze(Output.FCellCorrected(1:toneonsetframe(i),i,:)));
        repbaseline = repmat(baseline,[size(Output.FCellCorrected,1) 1]);
        DeltaFCellCorrected(:,i,:) = 100*(squeeze(Output.FCellCorrected(:,i,:))-repbaseline)./repbaseline;
        baseline = nanmean(squeeze(Output.FCellCorrectedSmooth(1:toneonsetframe(i),i,:)));
        repbaseline = repmat(baseline,[size(Output.FCellCorrectedSmooth,1) 1]);
        DeltaFCellCorrectedSmooth(:,i,:) = 100*(squeeze(Output.FCellCorrectedSmooth(:,i,:))-repbaseline)./repbaseline;
    end
    Output.DeltaFCellCorrected = DeltaFCellCorrected;
    Output.brightcells = brightcells;
    Output.DeltaFCellCorrectedSmooth = DeltaFCellCorrectedSmooth;
    [x y z] = size(Output.DeltaFCellCorrectedSmooth);
    %Find deconvovled spikes
    F = Output.fluoAllCorr';
    [nrois, nfr] = size(F);
    %Inidicates exponential decay used in deconvolution. Plenz lab uses
    %1.5, here using default
    tau = [];%ones(1,nrois)*1.5;
    %firing rate = lam*dt. Plenz lab uses 1. Here, using default.
    lam = [];
    %sig is a metric of F deviation from mean (median average deviation).
    %the MAD was multiplied by the constant 1.4826, to approximate std
    %estimation for normally distributed data, but the data are not, so it
    %is removed.
    sig = (mad(F,1,2)');%*1.4826;
    t = linspace(T,T*nfr,nfr);
    spikes = deconSpikesFromCal(F,t,T,sig,lam,tau);
    spikes.tau = tau;
    spikes.dt = T;
    %Reorganize spike trains into rasters
    %Spike detection threshold
    dconthresh = 3*mad(spikes.vals);
    %get raster from spikes
    Output.SpikeRaster = getROIraster(spikes,dconthresh,nrois,nfr,t,T);
    S=nan(x,y,z);
    Sthresh=nan(x,y,z);
    %trial taper to deal with edge effects from non-continuous DAQ
    win = hann(10)';
    Output.SpikeRaster.winlen = length(win)/2;
    win = [win(1:length(win)/2) ones(1,x-length(win)) win(1+length(win)/2:end)]; %make sure window is equal to length of TimingInfo.SeqEndVals
    for iii = 1:y
        for ii = 1:z
            try 
                S(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
                win.*Output.SpikeRaster.S(ii,TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2));
            catch
                temp_win = win
                temp_win(end-Output.SpikeRaster.winlen)=[]
                S(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
                temp_win.*Output.SpikeRaster.S(ii,TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2));
            end
            try
            Sthresh(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
            win.*Output.SpikeRaster.Sthresh(ii,TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2));
            catch
                temp_win = win
                temp_win(end-Output.SpikeRaster.winlen)=[]
                S(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
                temp_win.*Output.SpikeRaster.S(ii,TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2));
            end
        end
    end
    Output.SpikeRaster.S = S;
    Output.SpikeRaster.Sthresh = Sthresh;
    %Find S
    S=[];
    Sthresh=[];
    for i = 1:size(Freqs,2)
        baseline = nanmean(squeeze(Output.SpikeRaster.S(1:toneonsetframe(i),i,:)));
        repbaseline = repmat(baseline,[size(Output.SpikeRaster.S,1) 1]);
        S(:,i,:) = 100*(squeeze(Output.SpikeRaster.S(:,i,:))-repbaseline)./repbaseline;
        baseline = nanmean(squeeze(Output.SpikeRaster.Sthresh(1:toneonsetframe(i),i,:)));
        repbaseline = repmat(baseline,[size(Output.SpikeRaster.Sthresh,1) 1]);
        Sthresh(:,i,:) = 100*(squeeze(Output.SpikeRaster.Sthresh(:,i,:))-repbaseline)./repbaseline;
    end
    Output.SpikeRaster.DeltaS = S;
    Output.SpikeRaster.DeltaSthresh = Sthresh;
    save([LoadPath '/Fluorescence.mat'], 'Output')
end
