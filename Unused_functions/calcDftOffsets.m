function [offsets imTemplate] = calcDftOffsets(IMG,dftResolution,imTemplate)
% IMG should be MxNxFrames double

% Makes use of "dftregistration" code package from:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

%% FIND SEGMENT OF MOVIE WITH HIGH CORR VALUES
winsize = 99;
nframes = size(IMG,3);
sampInd = round(linspace(1,nframes,50));
corrSeq = zeros(length(sampInd),3);

tstartCorr = tic;

if ~exist('imTemplate','var')
    disp('Finding "stable" segment of movie....');
    for k = 1:length(sampInd)
        lwin = max(1, sampInd(k)-winsize);
        rwin = min(sampInd(k)+winsize, nframes);
        imgSeqTmp = reshape(IMG(:,:,lwin:rwin), size(IMG,1)*size(IMG,2),size(IMG(:,:,lwin:rwin),3));
        rho = corr(imgSeqTmp);
        corrSeq(k,:) = [mean(rho(:)) lwin rwin];
    end %K
    clear imgSeqTmp
    
    pkCorrInd = find(corrSeq(1:end-1,1) == max(corrSeq(1:end-1,1)));
    tEndCorr = toc(tstartCorr);
    disp([' Time to find stable portion of movie ' ,num2str(tEndCorr)/60,' minutes'])
    
    
   
    timPtWindow = (corrSeq(pkCorrInd,3) + corrSeq(pkCorrInd,2))/2;
    I = (mean(IMG(:,:,timPtWindow-50:timPtWindow+50),3));
    fixed = (I - min(I(:)))./range(I(:));
    imTemplate = fft2(fixed);
end
if isempty(gcp('nocreate'))
    parpool('local');
end
 %% perform DFT registration & get motion correction offsets
disp('')
disp('Finding motion correction coordinates....');
tStartMotionOffsets = tic;
parfor j = 1 :  nframes
    % using Fourier transformation of images for registration
    error = dftregistration(imTemplate,fft2(IMG(:,:,j)),dftResolution);
    ty(j) = error(3);
    tx(j) = error(4);
end
offsets = [ty' tx'];
telapsedMotionOffset = toc(tStartMotionOffsets);
disp('')
disp(['     Time Elapsed for Motion Offsets was: ', num2str(telapsedMotionOffset), ' seconds'])

end