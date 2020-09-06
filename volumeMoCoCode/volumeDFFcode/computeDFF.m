%% Compute dF/F Traces from raw images and defined cell centers
% Algorithms for ROI identification and Fluo trace extraction written by
% Dan Winkowski (DW)
% Zac Bowen put together flexible function and did all further edits
function [movID] = computeDFF(expDir,exptName,RegFile,exptVars,winSizeSeconds,percentBaselineSub,neuropilSubPercent)

narginchk(4,7)
if ~exist('winSizeSeconds','var'); winSizeSeconds = 10; end % default to 10 second window
if ~exist('percentBaselineSub','var'); percentBaselineSub = 50; end % default to 50% baseline subtraction
if ~exist('neuropilSubPercent','var'); neuropilSubPercent = 70; end % default to 70% neuropil subtraction

total_start = tic;
%% Dan experimental parameters
cd(expDir)

movID = ['DFF_' exptName];
datafile = [fullfile(expDir,movID) '.mat'];
if ~exist(datafile,'file')
    disp(['Creating data file', datafile]);
    save(datafile,'movID','-v7.3');
end

warning('off','stats:pvaluedw:ExactUnavailable')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DW05062015 - added to make the ROI ring more flexible and dependent on mag
%size rather than fixed
expectedNeuronDiamMicrons = 10; %expected diameter of neuron in microns
expectedNeuropilDiamMicrons = 20; %expected diameter of neuropil  in microns
expectedNeuronRadiusPix = round(expectedNeuronDiamMicrons/exptVars.micronsPerPixel)/2;
expectedNeuropilRadiusPix = round(expectedNeuropilDiamMicrons/exptVars.micronsPerPixel)/2; %need to remove outer boundary of neurons
% percNP = 0.9; %DW 11232015- added for Npil correction - percentage of neuropil signal to substract (Svoboda uses 1; Komiyama uses 0.9; Kerlin/Reid use 0.7 - hmmmmm......)
fractionNP = neuropilSubPercent / 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read RAW file into workspace
fullpathIMG = fullfile(expDir,RegFile);
fh = fopen(fullpathIMG);
% imgfile = fread(fh , inf , 'uint16=>uint16');
imgfile = fread(fh , exptVars.numVolumes*exptVars.segmentSize/2 , 'uint16=>uint16');
IMG = reshape(imgfile,[exptVars.dimX, exptVars.dimY, exptVars.numVolumes]);
% IMG = permute(IMG, [2 1 3]); % already reshaped in the motion correction code
fclose(fh);
clear imgfile
meanIMG = mean(IMG,3); % Mean image for cell center clicking
norm_meanIMG = (meanIMG - repmat(min(meanIMG(:)),size(meanIMG))) ./ repmat(range(meanIMG(:)),size(meanIMG));

save(datafile,'meanIMG','norm_meanIMG','-append');

%% Select neuron centers - FRAUNHOFER -- May have to make xc/yc an input to computeDFF given the Platform layout
% function [xc,yc] = getCellCenters(exptDir,exptName,meanIMG?)
cDefFile = dir(['CellDefFile_' exptName '.mat']);
if ~isempty(cDefFile)
%     xyCoordsTmp = load([expDir ,'/',cDefFile.name]);
% %     xc = xyCoordsTmp.CellCenters(:,1);
% %     yc = xyCoordsTmp.CellCenters(:,2);
%     yc = xyCoordsTmp.ptsIdx(:,2);
%     xc = xyCoordsTmp.ptsIdx(:,3);
    load(['CellDefFile_' exptName '.mat']) %updated for volume analysis (since multiple cdefs in 1 dir
else
    figure; imagesc(meanIMG); colormap('gray'); axis('square')
    SelectText = 'Click on Neuron Centers, please.........';
    disp ( SelectText );
    [xc, yc] = getpts; %  manually select centers of the neurons
    CellCenters(:,1) = xc;
    CellCenters(:,2) = yc;
    save(['CellDefFile_' exptName '.mat'],'CellCenters');
    %     set(h1,'position',hPos1);
    return
end

% centsel_finish = toc(centsel_start);
% fprintf('Selecting centers took %.1f minutes\n',centsel_finish/60)


%% PREALLOCATE
traceExtract_start = tic;

pcimg = cell (length(xc) , 13 , 360 );
imgCrop = cell ( length(xc) , 27, 27 );
imgCropNorm = cell (length(xc), 28 , 28 );
roiBoundaries = cell ( length(xc) , 360 , 3 );
smRoiBoundaries = cell ( length(xc) , 360 , 3 );
ROIOut = zeros ( 360 , 2 );
ROIIn =  zeros ( 360 , 2 );
ROIxvOut = cell ( length(xc) , 360 , 1 );
ROIyvOut = cell ( length(xc) , 360 , 1 );
ROIxvIn = cell ( length(xc) , 360 , 1 );
ROIyvIn = cell ( length(xc) , 360 , 1 );
roiBWout = cell ( length(xc) , size(norm_meanIMG,1) , size(norm_meanIMG,2));
roiBWin = cell ( length(xc) , size(norm_meanIMG,1) , size(norm_meanIMG,2));
roiBW2 = cell ( length(xc) , size(norm_meanIMG,1) , size(norm_meanIMG,2));
%     xvOut = cell ( length(xc) , 360 , 1 );
%     yvOut = cell ( length(xc) , 360 , 1 );
%     xvIn = cell ( length(xc) , 360 , 1 );
%     yvIn = cell ( length(xc) , 360 , 1 );
%     BWout = cell ( length(xc) , size(Imc_norm,1) , size(Imc_norm,2));
%     BWin = cell ( length(xc) , size(Imc_norm,1) , size(Imc_norm,2));
%     BW2 = cell ( length(xc) , size(Imc_norm,1) , size(Imc_norm,2));

%neuropil
npBoundaries = cell ( length(xc) , 360 , 3 );
smNpBoundaries = cell ( length(xc) , 360 , 3 );
NeuropilOut = zeros ( 360 , 2 );
NeuropilIn =  zeros ( 360 , 2 );
NPxvOut = cell ( length(xc) , 360 , 1 );
NPyvOut = cell ( length(xc) , 360 , 1 );
NPxvIn = cell ( length(xc) , 360 , 1 );
NPyvIn = cell ( length(xc) , 360 , 1 );
npBWout = cell ( length(xc) , size(meanIMG,1) , size(meanIMG,2));
npBWin = cell ( length(xc) , size(meanIMG,1) , size(meanIMG,2));
npBW2 = cell ( length(xc) , size(meanIMG,1) , size(meanIMG,2));

%% PREALLOCATE FLUO AND NPFLUO MATRICES

rawFluo = zeros( size(IMG,3) , length(xc) );
npFluo = zeros( size(IMG,3) , length(xc) );
npSubFluoSmooth = zeros( size(IMG,3) , length(xc) );

% parsedFluoTraces = zeros( niters , tarFnums+1 , size(xc,1) );
%parsedDFF =  cell(1 , nstims(1)*nstims(2) )  ;% , tarFnums+1 , size(xc,1) );

%% FIND THE BOUNDARIES OF CLICKED NEURONS
for pp = 1:length(xc)
    xpt=xc(pp);
    ypt=yc(pp);
    imgCrop{pp} = imcrop(meanIMG,[xpt-15 ypt-15 31 31]);
    imgCropNorm{pp} = (imgCrop{pp} - min(imgCrop{pp}(:)))  ./ (max(imgCrop{pp}(:)) - min(imgCrop{pp}(:)));
    pcimg{pp} = imgpolarcoord (imgCropNorm{pp} );  % this comes from Matlab Central Download
    
    RingPks = zeros(size(pcimg{pp},2),1); %reset vals to zero
    RingNeuropil = zeros(size(pcimg{pp},2),1); %reset vals to zero
    
    tmpNeuron = pcimg{pp}; %iterate this for each selected neuron
    
    for cc = 1:size(tmpNeuron,2) % for every direction - find the inner part of the ring
        
        pkTmp = find(diff(tmpNeuron(:,cc)) == min(diff(tmpNeuron(:,cc)))); % DW07122015_changed to make this more robust - seems to be working right now - continue testing
        
        if ~isempty(pkTmp)
            if length(pkTmp) > 1 %more than one pixel identified - grab the first one
                if pkTmp(1) < expectedNeuronRadiusPix && pkTmp(1) > 2
                    RingPks(cc) = pkTmp(1);
                else
                    RingPks(cc) = expectedNeuronRadiusPix;
                end
            else
                if pkTmp < expectedNeuronRadiusPix && pkTmp(1) > 2
                    RingPks(cc) = pkTmp;
                else
                    RingPks(cc) = expectedNeuronRadiusPix;
                end
            end
        elseif cc == 1 % if its'the first direction and no peaks are found
            RingPks(cc) = expectedNeuronRadiusPix;  %made this dependent on mag factor DW_02022015
        else
            RingPks(cc) = RingPks(cc-1);
        end
        ROIOut(cc,:) = [ deg2rad(cc)  RingPks(cc) ];
        ROIIn(cc,:) = [ deg2rad(cc)   RingPks(cc)-2 ];
        
        %%%%%%NEED TO INCLUDE NEUROPIL SIGNAL (SVOBODA LAB USES ~20 UM
        %%%%%%FROM CELL CTR EXCLUDING ROIs ALL OTHER NEURONS - NATURE
        %%%%%%2015)
        NeuropilIn(cc,:) = [deg2rad(cc) ROIOut(cc,2)+1];
        NeuropilOut(cc,:) = [deg2rad(cc) expectedNeuropilRadiusPix];
        
        %%%%%%NEED TO INCLUDE NEUROPIL SIGNAL (SVOBODA LAB USES ~20 UM
        %%%%%%FROM CELL CTR EXCLUDING ROIs ALL OTHER NEURONS - NATURE
        %%%%%%2015)
    end
    
    roiBoundaries{pp} = [ ROIIn(:,1) ROIIn(:,2)  ROIOut(:,2) ]; % [PolarCoords (0-2Pi)    InnerRing     OuterRing]
    smRoiBoundaries{pp} = [ ROIIn(:,1) smooth(ROIIn(:,2),10)  smooth(ROIOut(:,2),10) ]; % [PolarCoords (0-2pi)    InnerRing     OuterRing]
    %DW 11232015 - included neuropil variables
    npBoundaries{pp} = [ NeuropilIn(:,1) NeuropilIn(:,2) NeuropilOut(:,2) ];
    smNpBoundaries{pp} = [ NeuropilIn(:,1) smooth(NeuropilIn(:,2),10) smooth(NeuropilOut(:,2),10) ];
    
    % CREATE MASKS FOR ALL CLICKED ROIS, THEN SHOW THEM -- DW 11232015
    % renamed variable for consistency
    ROIxvOut{pp} =  xpt + smRoiBoundaries{pp}(:,3) .* (cos(smRoiBoundaries{pp}(:,1))) ;
    ROIyvOut{pp} =  ypt + smRoiBoundaries{pp}(:,3) .* (sin(smRoiBoundaries{pp}(:,1))) ;
    ROIxvIn{pp} =  xpt + smRoiBoundaries{pp}(:,2) .* (cos(smRoiBoundaries{pp}(:,1))) ;
    ROIyvIn{pp} =  ypt + smRoiBoundaries{pp}(:,2) .* (sin(smRoiBoundaries{pp}(:,1))) ;
    roiBWout{pp} = poly2mask( ROIxvOut{pp} , ROIyvOut{pp} , size(meanIMG,1) , size(meanIMG,2));
    roiBWin{pp} = poly2mask( ROIxvIn{pp} , ROIyvIn{pp} , size(meanIMG,1) , size(meanIMG,2));
    roiBW2{pp} =  roiBWout{pp} -  roiBWin{pp};
    if sum(roiBW2{pp}(:) < 0) > 0 %accounts for inner diameter extending beyond outer diameter
        roiBW2 {pp} (roiBW2 {pp} < 0 ) = 0;
    end
    
    %DW 11232015 - included for neuropil correction
    % CREATE MASKS FOR NEUROPIL CLICKED ROIS, THEN SHOW THEM
    NPxvOut{pp} =  xpt + smNpBoundaries{pp}(:,3) .* (cos(smNpBoundaries{pp}(:,1))) ;
    NPyvOut{pp} =  ypt + smNpBoundaries{pp}(:,3) .* (sin(smNpBoundaries{pp}(:,1))) ;
    NPxvIn{pp} =  xpt + smNpBoundaries{pp}(:,2) .* (cos(smNpBoundaries{pp}(:,1))) ;
    NPyvIn{pp} =  ypt + smNpBoundaries{pp}(:,2) .* (sin(smNpBoundaries{pp}(:,1))) ;
    npBWout{pp} = poly2mask( NPxvOut{pp} , NPyvOut{pp} , size(meanIMG,1) , size(meanIMG,2));
    npBWin{pp} = poly2mask( NPxvIn{pp} , NPyvIn{pp} , size(meanIMG,1) , size(meanIMG,2));
    npBW2{pp} =  npBWout{pp} -  npBWin{pp};
    if sum(npBW2{pp}(:) < 0) > 0 %accounts for inner diameter extending beyond outer diameter
        npBW2 {pp} (npBW2 {pp} < 0 ) = 0;
    end
end

%DW 11232015 - adjusted for inclusion of neuropil correction
% correct for overlapping ROIs (exclude from both)
disp('Adjusting ROI masks for overlap....');
tStartROICorr = tic;
AllMasksTMP =  sum ( cat ( 3 , roiBWout{:} ) , 3 ); % first term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
[oLapRoiY, oLapRoiX] = find( AllMasksTMP > 1 );
for ii = 1:pp
    for yy = 1:length(oLapRoiX)
        roiBWout{ii}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
        roiBW2{ii}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
    end
end
ROImap =  sum ( cat ( 3 , roiBW2{:} ) , 3 ); % first term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
tElapsedROICorr = toc(tStartROICorr);
disp(['    Time elapsed for ROI mask Correction was ', num2str(tElapsedROICorr/60),' minutes']);

% correct for neuropil overlap with ROIs and
disp('Adjusting neuropil masks for overlap....');
tStartNPilCorr = tic;
for nn = 1:size(npBWout,1)
    tmpNp=[];
    tmpNp = npBWout{nn};
    for rr = 1:size(roiBWout,1)
        tmp = roiBWout{rr} + tmpNp;
        [oLapNPY,oLapNPX] = find( tmp > 1 );
        for yy = 1:length(oLapNPX)
            npBWout{nn}(oLapNPY(yy),oLapNPX(yy)) = 0;
            npBWin{nn}(oLapNPY(yy),oLapNPX(yy)) = 0;
        end
        for yy = 1:length(oLapRoiX)
            npBWout{nn}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
            npBWin{nn}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
        end
        
    end
end
AllNpilMasks =  sum ( cat ( 3 , npBWout{:} ) , 3 ); % first term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
tElapsedNPilCorr = toc(tStartNPilCorr);
disp(['    Time elapsed for Neuropil Mask Correction was ', num2str(tElapsedNPilCorr/60),' minutes']);

RoiNpMask(:,:,1) = zeros(size(AllNpilMasks));
RoiNpMask(:,:,2) = ROImap ;
RoiNpMask(:,:,3) = AllNpilMasks;

save(datafile, 'ROImap', 'xc','yc','-append');

traceExtract_finish = toc(traceExtract_start);
fprintf('Trace extraction took %.1f minutes\n',traceExtract_finish/60)

%% CALCULATE FLUORESCENCE TRACES for ALL STIMS

calcF_start = tic;

for nn = 1:length(xc)
    [r,c]=find(roiBW2{nn}~=0);
    
    tmp = NaN(length(r),size(IMG,3));
    for i = 1:length(r)
        tmp(i,:) = IMG(r(i),c(i),:);
    end
    rawFluo(:,nn) = nanmean(tmp,1);
    
end
clear roiBW2
calcF_finish = toc(calcF_start);
fprintf('Calculating fluorescence traces took %.1f minutes\n',calcF_finish/60)

%DW 11232015 - included Neuropil calculation
disp('Calculating Fluo traces for neuropil.....');
tStartfcalcNP = tic;
%DW11132015 - apply to neuropil calculation.....
for nn = 1:length(xc)
    [rNp,cNp] = find(npBWout{nn} ~= 0 );
    
    tmp = NaN(length(rNp),size(IMG,3));
    for i = 1:length(rNp)
        tmp(i,:) = IMG(rNp(i),cNp(i),:);
    end
    npFluo(:,nn) = nanmean(tmp,1);
end
tElapsedFcalcNP = toc(tStartfcalcNP);
disp(['    Time Elapsed for calculating fluo for neuropil was ',num2str(tElapsedFcalcNP/60),' minutes']);

clear IMG npBWout

npSubFluo = rawFluo - (fractionNP * npFluo);
% fluoAllCorr = fluoAllRaw - (percNP * NPfluoAll);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PUT THE DE-NOISING, PARSING AND DFF CALCULATIONS HERE - %% DW05192016 added for trialBased analyses

%% apply assymetrical (exponential) filter to traces - smoothing - MODIFIED FROM PLENZ LAB AVALANCHE CODE
[nrois , ~] = size(npSubFluo');
T = 1/exptVars.frameRate;

% win = 10; %Change based on fast rise and exponential decay
% tau = 1.5;
% wf = 0:win;
% wf = exp(-tau*T*wf);

% Normalizing raw F to get rid of negative vals before dF/F calc % Jul2020
adjF = npSubFluo - min(npSubFluo(:));
normF = adjF ./ max(adjF(:));

%smoothed raw dR/R or dF/F - commented out 9/16, uncommented Jun2020 for slower volume data
if T >= 0.2
    smoothwin = 3;
elseif T >= 0.13 && T <= 0.2
    smoothwin = 5; %7
elseif T <=0.13
    smoothwin = 9; %13
end
symmFLAG = 1; % 1 ->symmetric, 2->asym, T = frame time, Tau -> 1.5
% smoothwin = 9; % Because T == 0.033

for j = 1 : nrois % NEED FOR IMPROVEMENT - ELIMINATE FOR LOOP?
    
    %%%%% SYMMETRIC (new as of Sept2016)
    tmp = normF(:,j)';
    npSubFluoSmooth(:,j) = smooth2(tmp,symmFLAG,smoothwin,T,1.5); %Maybe smooth after DFF calc?
    
end

% performing moving average baseline subtraction for DFF
winsize = round(exptVars.frameRate * winSizeSeconds);
percent = percentBaselineSub;
DFF.smooth = slideWinSub(npSubFluoSmooth',winsize,percent);
DFF.raw = slideWinSub(normF',winsize,percent);

%% clean up raw (remove bad ROIs)

% Streamlined the "cleaning" stage June2020
% Jul2020 update: might not be a necessary step anymore with normalization changes made above, saving regardless
[DFF.smoothClean,keptRois.smoothClean] = cleanTrace(DFF.smooth);
[DFF.rawClean,keptRois.rawClean] = cleanTrace(DFF.raw);

dffParams.winSizeSeconds = winSizeSeconds;
dffParams.percentBaselineSub = percentBaselineSub;
dffParams.neuropilSubPercent = neuropilSubPercent;

save(datafile,'npSubFluo','DFF','keptRois','dffParams','-append')

total_finish = toc(total_start);
fprintf('computeDFF.m completed successfully in %.1f minutes\n',total_finish/60);


end