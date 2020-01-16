function output = ExtractFluorescence_datalist(input)

savepath = input.path;

%ExtractFluorescence is a program that load cell coordinates and extract brightness rings for cell
%membranes and neuropil
%Nikolas Francis 2016
tic();
%Expected frame rate for imaging.
fps = input.expectedFPS;
bb=strsplit(input.path,'/'); % KA- changed strsep to strsplit
%Load cell definitions
load([input.path '/' 'CellDefinitions.mat'])
xc = ptsIdx(:,2);
yc = ptsIdx(:,3);
border = input.border:dimX;
bordermask = ones(dimX,dimX);
bordermask(border,border)=0;
for pp = 1:length(xc)
    xpt = xc(pp);
    ypt= yc(pp);
    %Cell ROIs
    ROIxvOut{pp} =  xpt + smRoiBoundaries{pp}(:,3) .* (cos(smRoiBoundaries{pp}(:,1))) ;
    ROIyvOut{pp} =  ypt + smRoiBoundaries{pp}(:,3) .* (sin(smRoiBoundaries{pp}(:,1))) ;
    ROIxvIn{pp} =  xpt + smRoiBoundaries{pp}(:,2) .* (cos(smRoiBoundaries{pp}(:,1))) ;
    ROIyvIn{pp} =  ypt + smRoiBoundaries{pp}(:,2) .* (sin(smRoiBoundaries{pp}(:,1))) ;
    roiBWout{pp} = poly2mask( ROIxvOut{pp} , ROIyvOut{pp} , dimX , dimX);
    roiBWin{pp} = poly2mask( ROIxvIn{pp} , ROIyvIn{pp} , dimX , dimX);
    roiBW2{pp} =  roiBWout{pp} -  roiBWin{pp};
    roiTOTAL{pp} =  roiBWout{pp} +  roiBWin{pp};
    %Account for inner diameter extending beyond outer diameter
    if sum(roiBW2{pp}(:) < 0) > 0
        warndlg('Inner and Outer ROIs overlapping!')
        roiBW2 {pp} (roiBW2 {pp} < 0 ) = 0;
    end
    %Neuropil ROIs
    NPxvOut{pp} =  xpt + smNpBoundaries{pp}(:,3) .* (cos(smNpBoundaries{pp}(:,1))) ;
    NPyvOut{pp} =  ypt + smNpBoundaries{pp}(:,3) .* (sin(smNpBoundaries{pp}(:,1))) ;
    NPxvIn{pp} =  xpt + smNpBoundaries{pp}(:,2) .* (cos(smNpBoundaries{pp}(:,1))) ;
    NPyvIn{pp} =  ypt + smNpBoundaries{pp}(:,2) .* (sin(smNpBoundaries{pp}(:,1))) ;
    npBWout{pp} = poly2mask( NPxvOut{pp} , NPyvOut{pp} , dimX , dimX);
    npBWin{pp} = poly2mask( NPxvIn{pp} , NPyvIn{pp} , dimX , dimX);
    npBW2{pp} =  npBWout{pp} -  npBWin{pp};
    %Account for inner diameter extending beyond outer diameter
    if sum(npBW2{pp}(:) < 0) > 0
        warndlg('Inner and Outer ROIs overlapping!')
        npBW2 {pp} (npBW2 {pp} < 0 ) = 0;
    end
end
%Correct for neuropil overlap with ROIs
disp('Adjusting NEUROPIL masks for overlap....');
AllMasksTMP = sum ( cat ( 3 , npBWin{:} ) , 3 );
AllMasksTMP = bordermask + AllMasksTMP + sum ( cat ( 3 , roiTOTAL{:} ) , 3 );
[oLapRoiY, oLapRoiX] = find( AllMasksTMP > 1 );
for ii = 1:pp
    for yy = 1:length(oLapRoiX)
        npBW2{ii}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
    end
end
%Correct for overlapping ROIs (exclude from both by setting values to 0)
disp('Adjusting ROI masks for overlap....');
%First term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
AllMasksTMP =  bordermask + sum ( cat ( 3 , roiBWout{:} ) , 3 );
[oLapRoiY, oLapRoiX] = find( AllMasksTMP > 1 );
for ii = 1:pp
    for yy = 1:length(oLapRoiX)
        roiBW2{ii}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
    end
end
%Load image sequences
opts = get_options_from_xml([input.path '/Experiment.xml']);
opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};
try
    GreenChannel = memmapfile([input.path '/shifted_img.raw'], ...
    'Format', opts.format, 'Repeat', opts.numframes);
catch
    GreenChannel = memmapfile([input.path '/greenchannelregistered.raw'], ...
    'Format', opts.format, 'Repeat', opts.numframes);
    disp(['check ' input.path 'for right raw file']);
end
%Find fluorescence for each ROI and neuropil. fluoAllRaw: ROI
%uncorrected for neuropil; NPfluoAll: Neuropil; fluoAllCorr:
%corrected ROI;
%Border Mask
fluoAllRaw(1:opts.numframes,length(xc))=0;
NPfluoAll(1:opts.numframes,length(xc))=0;
NPfluoAll_100pct(1:opts.numframes,length(xc))=0;
%Check if need to process movie in chunks or at once
if opts.numframes < input.maxframechunk
    greenChanImg=[];
    for frame = 1:opts.numframes
        greenChanImg(:,:,frame) = GreenChannel.Data(frame).channels;
        %fprintf('loading frame %d/%d/n', frame, opts.numframes);
    end
    disp('Calculating fluorescence traces for SOMA and NEUROPIL.....');
    for nn = 1:length(xc)
        fprintf('Processing cell %d/%d/n', nn, length(xc));
        [ r , c ] = find(roiBW2{nn} ~= 0 );
        if ~isempty(r)
            [ rNp , cNp ] = find(npBW2{nn} ~= 0 );
            for frame = 1:opts.numframes
                g = squeeze(greenChanImg(:,:,frame));
                fluoAllRaw(frame,nn) = nanmean(g(sub2ind(size(g),r,c)));
                NP = g(sub2ind(size(g),rNp,cNp));
                prtct = prctile(NP,80);
                NP = NP(NP<prtct);
                NPfluoAll(frame) = nanmean(NP);
            end
        end
    end
    fluoAllCorr = fluoAllRaw - (input.percNP * NPfluoAll);
else
    frameidx = [0 floor(opts.numframes/2) opts.numframes];
    for f = 1:length(frameidx)-1
        fidx = frameidx(f)+1:frameidx(f+1);
        greenChanImg=[];
        for frame = 1:length(fidx)
            greenChanImg(:,:,frame) = GreenChannel.Data(fidx(frame)).channels;
            %fprintf('loading frame %d/%d/n', fidx(frame), opts.numframes);
        end
        %Find flouro traces
        disp('Calculating fluorescence traces for SOMA and NEUROPIL.....');
        for nn = 1:length(xc)
            fprintf('Processing cell %d/%d/n', nn, length(xc));
            [ r , c ] = find(roiBW2{nn} ~= 0 );
            if ~isempty(r)
                [ rNp , cNp ] = find(npBW2{nn} ~= 0 );
                fluotemp=zeros(length(fidx),1);
                NPfluotemp=zeros(length(fidx),1);
                
                for frame = 1:length(fidx)
                    g = squeeze(greenChanImg(:,:,frame));
                    fluotemp(frame) = nanmean(g(sub2ind(size(g),r,c)));
                    NP = g(sub2ind(size(g),rNp,cNp));
                    prtct = prctile(NP,80);
                    NP = NP(NP<prtct);
                    if any(isnan(NP)) == 1
                        NP = 0
                        disp(['caught a NaN in ' input.path '-check it out' fidx(frame)])
                    else
                    end
                    NPfluotemp(frame) = nanmean(NP);
                    NPfluotemp_100pct(frame) = nanmean(g(sub2ind(size(g),rNp,cNp)));
                end
                fluoAllRaw(fidx,nn) = fluotemp;
                NPfluoAll(fidx,nn) = NPfluotemp;
                NPfluoAll_100pct(fidx,nn) = NPfluotemp_100pct;
            end
        end
        fluoAllCorr = fluoAllRaw - (input.percNP * NPfluoAll);
        fluoAllCorr_100pct = fluoAllRaw - (input.percNP * NPfluoAll_100pct);
    end
end
%Extract timing parameters
PsigFilename = dir([input.path,'/ART*.mat'])

PsignalFile=[input.path '/' PsigFilename.name];
ThorFile = [input.path '/timing.txt'];
TimingInfo = getTimingInfo(ThorFile, PsignalFile, fps);
%Parse flourscnence into sample X trial X cell. sample values are
%initialized with nan, and # samples is determined by the expected
%trial length from TimingInfo. If trial data were shorter than the
%expected amount, then the remaining samples will be nan.
FCell=nan(TimingInfo.tarFnums,length(TimingInfo.SeqEndVals),size(fluoAllRaw,2));
FNeuropil=nan(TimingInfo.tarFnums,length(TimingInfo.SeqEndVals),size(fluoAllRaw,2));
FCellCorrected=nan(TimingInfo.tarFnums,length(TimingInfo.SeqEndVals),size(fluoAllRaw,2));
for iii = 1:length(TimingInfo.SeqEndVals)
    for ii = 1:size(fluoAllRaw,2)
            try
            FCell(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
            fluoAllRaw(TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2),ii);
            FNeuropil(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
            NPfluoAll(TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2),ii);
            FCellCorrected(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
            fluoAllCorr(TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2),ii);
            catch
            frame_check = TimingInfo.FrameIdx(iii,2)-TimingInfo.FrameIdx(iii,1)+1;    
            FCell(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii) frame_check ]),iii,ii) = ...
            fluoAllRaw(TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2),ii);
            FNeuropil(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii) frame_check]),iii,ii) = ...
            NPfluoAll(TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2),ii);
            FCellCorrected(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii) frame_check]),iii,ii) = ...
            fluoAllCorr(TimingInfo.FrameIdx(iii,1):TimingInfo.FrameIdx(iii,2),ii);
            end
    end
end
%Save traces
Output.FCell = FCell;
Output.FNeuropil = FNeuropil;
Output.FCellCorrected = FCellCorrected;
Output.fluoAllRaw  = fluoAllRaw;
Output.NPfluoAll = NPfluoAll;
Output.NPfluoAll_100pct= NPfluoAll_100pct;
Output.fluoAllCorr = fluoAllCorr;
Output.fluoAllCorr_100pct = fluoAllCorr_100pct;
Output.TimingInfo = TimingInfo;
save([input.path '/Fluorescence.mat'], 'Output')
%Plot brightness stats and average segemented cell
if input.plotbrightness
    figure
    subplot(2,2,1:2)
    %Find if cell brightness is > 0 and plot each cells brightness, ie.
    %cellular ROI vs. Neuropil ROI (%)
    cellbrightness=squeeze(100*(nanmean(nanmean(FCell,2),1)- ...
        nanmean(nanmean(FNeuropil,2),1))./nanmean(nanmean(FNeuropil,2),1));
    for c = 1:length(cellbrightness)
        hold on
        if cellbrightness(c) >= 0
            bar(c,cellbrightness(c),'r')
        elseif cellbrightness(c) < 0
            bar(c,cellbrightness(c),'b')
        end
    end
    hold on
    aa=axis;
    text(2,aa(4)-1,['nanmean = ' num2str(nanmean(cellbrightness(cellbrightness>=0))) ...
        ';N=' num2str(sum(cellbrightness>=0)) ...
        ';3% N=' num2str(sum(cellbrightness>=3))],'color','r')
    plot([aa(1) aa(2)],[nanmean(cellbrightness(cellbrightness>=0)) ...
        nanmean(cellbrightness(cellbrightness>=0))],'r')
    text(2,aa(4)-5,['nanmean = ' num2str(nanmean(cellbrightness(cellbrightness<0))) ...
        ';N=' num2str(sum(cellbrightness<0))],'color','b')
    plot([aa(1) aa(2)],[nanmean(cellbrightness(cellbrightness<0)) ...
        nanmean(cellbrightness(cellbrightness<0))],'b')
    plot([aa(1) aa(2)],[3 3],'k--')
    ylabel([{'Cell brightness'};{'(% re. neuropil background)'}])
    title([input.path],'Interpreter','none')
    xlabel('Cell #')
    muIMG = squeeze(nanmean(greenChanImg,3))';
    AllPosCells=0;
    AllNegCells=0;
    ccount=0;
    %Plot average cell for bright and dim cells
    input.cellcropdim = 80;
    for pp = 1:length(xc)
        %Crop the neuron
        xpt=xc(pp);
        ypt=yc(pp);
        imgCrop = imcrop(muIMG,[xpt-(input.cellcropdim/2) ypt-(input.cellcropdim/2) ...
            input.cellcropdim input.cellcropdim]);
        if size(imgCrop,1) == input.cellcropdim+1 && size(imgCrop,2) == input.cellcropdim+1
            ccount = ccount+1;
            imgCropNorm = (imgCrop-min(imgCrop(:)))./(range(imgCrop(:)));
            if cellbrightness(pp) >= 0
                AllPosCells = AllPosCells + imgCropNorm;
            elseif cellbrightness(pp) < 0
                AllNegCells = AllNegCells + imgCropNorm;
            end
        end
    end
    AllPosCells=AllPosCells./ccount;
    AllNegCells=AllNegCells./ccount;
    clim=[min([AllPosCells(:); AllNegCells(:)]) max([AllPosCells(:); AllNegCells(:)])];
    subplot(2,2,3)
    imshow(AllPosCells,[clim])
    title([{'Average of Bright Rings'};{'Ring>Neuropil'}])
    subplot(2,2,4)
    title('Average of Dim Rings')
    imshow(AllNegCells,[clim])
    title([{'Average of Dim Rings'};{'Ring<Neuropil'}])
end
fprintf('Elapsed time: %g. minutes/n', toc()/60);
end
