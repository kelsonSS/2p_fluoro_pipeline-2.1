function Output = ExtractFluorescence(input,debug)
%ExtractFluorescence is a program that load cell coordinates and extract brightness rings for cell
%membranes and neuropil
%Nikolas Francis 2016
tic();
%Expected frame rate for imaging.
Output.Errors = {};
fps = input.expectedFPS;
bb=strsplit(input.path,'\'); 
LocalPath = fullfile(input.savepath ,  bb{end}) ;
%Load cell definitions
for e = 1:length(input.expname)
    filecheck = char(fullfile(LocalPath, input.expname(e),...
                     'Fluorescence.mat'));
   if exist (filecheck,'file') && ( ~exist('debug','var')||debug == 0)
        fprintf('\n file %s already registered',filecheck);
       continue
    end
    
try
    load(fullfile(LocalPath , char(input.expname{e}) , 'CellDefinitions.mat'));
catch 
    warning(' could not find %s. Analyzing next file',fullfile(bb{end},input.expname{e}));  
    Output.Errors(end+1,1)={[bb{end} input.expname{e}]};
    continue
end

TimingFile = fullfile(LocalPath, input.expname{e} , 'TimingInfo.mat');
if exist(TimingFile,'file')
         
          load(TimingFile);    
       
      else 
          % hotfix- consider rework
         temp = input;
         temp.expname = temp.expname{e};
         TimingInfo =  ExtractTimingParams(input,1);
         clear temp

      end
     try
     tarFnum = TimingInfo.tarFnums;
   
     if TimingInfo.FrameIdx(1,1) == 0;
         TimingInfo.FrameIdx = TimingInfo.FrameIdx+1; % matlab uses 1-based indexing
     end 
     on_frame = TimingInfo.FrameIdx(:,1);
     off_frame = TimingInfo.FrameIdx(:,2);
   %  if estimated frames timings are off at all correct them 
     if  abs(mean(off_frame(2:end)-on_frame(2:end)+ 1) - tarFnum ) > 1e-14 
        
        for ii = 1: length(on_frame)
            %TimingFestimate = off_frame(ii)-on_frame(ii)+ 1;
            %frame_diff = TimingFestimate - tarFnum-1;
            off_frame(ii) =  on_frame(ii)+ TimingInfo.tarFnums - 1 ;
        end 
     end 
     catch 
         warning('TimingFile incorrect for %s. Analyzing next file',fullfile(bb{end},input.expname{e}));  
    Output.Errors(end+1,1)={[bb{end} input.expname{e}]};
         continue
     end 





 fprintf(' Extracting %s \n',LocalPath)
 fprintf('               %s \n',input.expname{e} );
try
    xc = ptsIdx(:,2);
    yc = ptsIdx(:,3);
    if length(xc) ~= length(smRoiBoundaries(:,1))
        xc = xc(1:length( smRoiBoundaries(:,1) ));
        yc = yc(1:length( smRoiBoundaries(:,1) ));
    end 
    border = input.border:dimX;
    bordermask = ones(dimX,dimX);
    bordermask(border,border)=0;
            
    for pp = 1:length(xc)
       % disp(pp)
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
    
     
    opts = get_options_from_xml(fullfile(LocalPath, input.expname{e}, 'Experiment.xml'));
    opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};
   
    %Find fluorescence for each ROI and neuropil. fluoAllRaw: ROI
    %uncorrected for neuropil; NPfluoAll: Neuropil; fluoAllCorr:
    %corrected ROI;
%    
%     %NPfluoAll_100pct(1:opts.numframes,length(xc))=0;
%     %Check if need to process movie in chunks or at once
%     if opts.numframes < input.maxframechunk
%      img_path = fopen(fullfile(LocalPath,input.expname{e},'greenchannelregistered.raw'));
%      greenChanImg =  fread(img_path,'uint16=>uint16');  
%      greenChanImg = reshape(greenChanImg,opts.dimX,opts.dimY,[]);
%      numframes = size(greenChanImg,3)
%      %Border Mask
%     fluoAllRaw(1:numframes,length(xc))=nan;
%     NPfluoAll(1:opts.numframes,length(xc))=nan;
%      
%     disp('Calculating fluorescence traces for SOMA and NEUROPIL.....');
%         for nn = 1:length(xc)
%             fprintf('Processing cell %d/%d \n', nn, length(xc));
%             [ r , c ] = find(roiBW2{nn} ~= 0 );
%             if ~isempty(r)
%                 [ rNp , cNp ] = find(npBW2{nn} ~= 0 );
%                 for frame = 1: size(greenChanImg,3)
%                   % disp(frame)
%                     g = squeeze(greenChanImg(:,:,frame));
%                     fluoAllRaw(frame,nn) = nanmean(g(sub2ind(size(g),r,c)));
%                     NP = g(sub2ind(size(g),rNp,cNp));
% %                      prtct = prctile(NP,80);
% %                      NP = NP(NP<prtct);
%                      NPfluoAll(frame) = nanmean(NP);
%                 end
%             end
%         end
%         data_idx = ~isnan(fluoAllRaw(:,1));
%         fluoAllRaw = fluoAllRaw(data_idx,:);
%         NPfluoAll = NPfluoAll(data_idx,:);
%         
%         fluoAllCorr = fluoAllRaw - (input.percNP * NPfluoAll);
%     else
        
fullpathIMG = fullfile(LocalPath, input.expname{e}, 'greenchannelregistered.raw');

fh = fopen(fullpathIMG); 
numFrames_total =  FindRawImgSize(fh,[opts.dimX opts.dimY]);
ItemsPerImage = opts.dimX* opts.dimY ;
chunkSize =  input.maxframechunk  ; % frames 
chunks = ceil(numFrames_total/chunkSize); 
fclose(fh);
NPfluoAll_temp = [];
fluoAllRaw_temp= [];

tstartCorr = tic;
for chunk_count = 1:chunks
         fh = fopen(fullpathIMG); 
         start_idx = (chunk_count-1) * ItemsPerImage * chunkSize * 2; %bytes  ;
    fseek(fh,start_idx,'bof');
    fprintf('%s : %d of %d \n', input.expname{1},chunk_count,chunks)
    % check to ensure we don't go over IMG size on last img
    if chunk_count == chunks 
        greenChanImg = fread(fh,inf,'uint16');
    else 
        greenChanImg =fread(fh,ItemsPerImage * chunkSize,'uint16');
    end 
    fclose(fh);
               
                try
                greenChanImg = reshape(greenChanImg,opts.dimX,opts.dimY,[]);
                catch 
                    % pad frame with zeros if it does not evenly divide for
                    % some reason 
                   divisor = opts.dimX*opts.dimY;
                   remainder = mod(length(greenChanImg),divisor);
                   greenChanImg = [greenChanImg;zeros( (divisor -remainder) ,1)];
                   greenChanImg = reshape(greenChanImg,opts.dimX,opts.dimY,[]);
                end 
                 greenChanImg = permute(greenChanImg,[2 1 3]);
                      
              numframes = size(greenChanImg,3);
            fidx = size(greenChanImg,3);
            %Border Mask
            fluoAllRaw(1:numframes,length(xc))=nan;
            NPfluoAll(1:numframes,length(xc))=nan;
           
          
            %Find flouro traces
            disp('Calculating fluorescence traces for SOMA and NEUROPIL.....');
            for nn = 1:length(xc)
                fprintf('Processing cell %d/%d \n', nn, length(xc));
                [ r , c ] = find(roiBW2{nn} ~= 0 );
                if ~isempty(r)
                    [ rNp , cNp ] = find(npBW2{nn} ~= 0 );
                   
                   % we have to directly index values for ND matricies so
                   % to do this nice and linear we will repeat our indicies
                   % by the number pixels in each frame this will get us to
                   % the same indicies in the next frame we then average
                   % those pixels to get the average fluorescence 
                   
                   imageInPixels = opts.dimX * opts.dimY;
                   cell_ind = sub2ind( [opts.dimX,opts.dimY],r,c);
                   numPixels = length(cell_ind);
                   cell_ind_full = [ cell_ind + imageInPixels *(0:numframes-1)];
                   g = greenChanImg(cell_ind_full);
                   fluoAllRaw_temp(:,nn) =  nanmean(g);
                  %repeat for neuropil
                   np_ind = sub2ind( [opts.dimX,opts.dimY],rNp,cNp);
                   numPixels = length(np_ind);
                   np_ind_full = [ np_ind + imageInPixels *(0:numframes-1)];
                   g = greenChanImg(np_ind_full);
                   NPfluoAll_temp(:,nn) =  nanmean(g);
                
                
                   
                   
                    
%                     Old reliable code
%                     fluotemp = nan(length(fidx),1);
%                     NPfluotemp = nan(length(fidx),1);
%                   
%                 
%                     
%                     for frame = 1:fidx
%                         g = squeeze(greenChanImg(:,:,frame));
%                         fluotemp(frame) = nanmean(g(sub2ind(size(g),r,c)));
%                         NP = g(sub2ind(size(g),rNp,cNp));
% %                         prtct = prctile(NP,80);
% %                         NP = NP(NP<prtct);
%                        if any(isnan(NP))
%                             NP = 0
%                         else
%                         end
%                         NPfluotemp(frame) = nanmean(NP);
%                         %NPfluotemp_100pct = mean(g(sub2ind(size(g),rNp,cNp)));
%                     end
%                     fluoAllRaw(:,nn) = fluotemp;
%                     NPfluoAll(:,nn) = NPfluotemp;
                    %%NPfluoAll_100pct(fidx,nn) = NPfluotemp_100pct;
                end
            end
  fluoAllRaw = cat(1,fluoAllRaw,fluoAllRaw_temp);
  NPfluoAll = cat(1,NPfluoAll,NPfluoAll_temp);
end
            fluoAllCorr = fluoAllRaw - (input.percNP * NPfluoAll);
            %fluoAllCorr_100pct = fluoAllRaw - (input.percNP * NPfluoAll_100pct);
    %Extract timing parameters
    
        try 
        copyfile(fullfile(input.path, input.expname{e}, input.psignalfiles{e}),...
                 fullfile(LocalPath,   input.expname{e}))
        catch 
            warning('could not move %s %s to new Analyzed folder,Skipping expt',...
                input.path,input.expname)
          continue
        end 
 
       
       
       
      
      
    %Parse flouresence into sample X trial X cell. sample values are
    %initialized with nan, and # samples is determined by the expected
    %trial length from TimingInfo. If trial data were shorter than the
    %expected amount, then the remaining samples will be nan.
  
   
     
     
    try 
    FCell=nan(TimingInfo.tarFnums,length(TimingInfo.SeqEndVals),size(fluoAllRaw,2));
    FNeuropil=nan(TimingInfo.tarFnums,length(TimingInfo.SeqEndVals),size(fluoAllRaw,2));
    FCellCorrected=nan(TimingInfo.tarFnums,length(TimingInfo.SeqEndVals),size(fluoAllRaw,2));
    for iii = 1:length(TimingInfo.SeqEndVals)
        for ii = 1:size(fluoAllRaw,2)
            if off_frame(iii) > size(fluoAllRaw,1)
                continue
            end 
      %     try
                %disp(iii)
            FCell(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
                fluoAllRaw(on_frame(iii):off_frame(iii),ii);
            FNeuropil(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
                NPfluoAll(on_frame(iii):off_frame(iii),ii);
            FCellCorrected(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii)]),iii,ii) = ...
                fluoAllCorr(on_frame(iii):off_frame(iii),ii);
            %  catch
               % if any(TimingInfo.FrameIdx(iii,:)> opts.numframes )
               %     continue
              %  end 
              
                       
              
%             frame_check = off_frame(iii)-on_frame(iii)+1;    
%             FCell(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii) frame_check ]),iii,ii) = ...
%             fluoAllRaw(on_frame(iii):off_frame(iii),ii);
%             FNeuropil(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii) frame_check]),iii,ii) = ...
%             NPfluoAll(on_frame(iii):off_frame(iii),ii);
%             FCellCorrected(1:min([TimingInfo.tarFnums TimingInfo.SeqEndVals(iii) frame_check]),iii,ii) = ...
%             fluoAllCorr(on_frame(iii):off_frame(iii),ii);
           %  end
        end
    end
    catch 
         warning('extraction incorrect for %s. Analyzing next file',fullfile(bb{end},input.expname{e}));  
    Output.Errors(end+1,1)={[bb{end} input.expname{e}]};
         continue
     end  
        
    
    %Save traces
    Output.FCell = FCell;
    Output.FNeuropil = FNeuropil;
    Output.FCellCorrected = FCellCorrected;
    Output.fluoAllRaw  = fluoAllRaw;
    Output.NPfluoAll = NPfluoAll;
    %Output.NPfluoAll_100pct= NPfluoAll_100pct;
    Output.fluoAllCorr = fluoAllCorr;
   % Output.fluoAllCorr_100pct = fluoAllCorr_100pct;
    Output.TimingInfo = TimingInfo;
    save(fullfile(LocalPath, input.expname{e}, 'Fluorescence.mat'), 'Output')
    %Plot brightness stats and average segemented cell
    if exist('debug','var')
       figure 
       ax= axes;
       plot(fluoAllCorr)
       for ii =1:length(on_frame)
           hold on
           plot([on_frame(ii) off_frame(ii) ], [ax.YLim(2) ax.YLim(2)])
       end 
       
       title(fullfile(input.path,input.expname{e}),'Interpreter','none')
        
       
       figure;
       imagesc(squeeze(nanmean(greenChanImg(:,:,1:300),3)));
       
       hold on 
       plot(xc,yc,'r.')
       
        
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
        title([input.path input.expname{e}],'Interpreter','none')
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

catch
       warning('Error: unable to extract Cell Fluorescence.Manual inspection needed.')
      warnpath = fopen(fullfile(input.savepath,bb{end},input.expname{e},'Warnings.txt'),...
                      'w+');
     fprintf(warnpath,'\n Error: unable to extract Cell Fluorescence.Manual inspection needed.  \n');
    
     fclose(warnpath); 
 end 
clear NPfluoAll
clear fluoAllRaw
clear fluoAllCorr
end
fprintf('Elapsed time: %g. minutes\n', toc()/60);

