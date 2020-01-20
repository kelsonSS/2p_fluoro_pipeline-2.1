function ExamineFluorescence(input,paths)
% this function quickly displays plots for for all fluorescence
% files that were extracted

for path_idx = 1:length(paths)
    
    in_path = paths{path_idx};
    bb=strsplit(in_path,filesep);
    file_path = fullfile(input.savepath ,  bb{end}) ;
    
    allfiles = dir([file_path , '\**\Fluorescence.*']);
    
    allfiles = struct2cell(allfiles);
    allfiles = allfiles(2,:)';
    
    if isempty(allfiles)
        continue
    end
    
    for file = 1:length(allfiles)
        try
            load(fullfile(allfiles{file},'Fluorescence.mat'))
            load(fullfile(allfiles{file},'TimingInfo.mat'))
        catch
            continue
        end
        FCellCorrected = Output.FCellCorrected;
        trialdur = size(FCellCorrected,1);
        trials   = size(FCellCorrected,2);
        neurons  = size(FCellCorrected,3);
        
        B_Vec = repmat(nanmean(FCellCorrected(1:30,:,:)),[trialdur,1,1]);
        Vec_DFF = (FCellCorrected -B_Vec)./B_Vec * 100;
        gausfilt = fspecial('gaussian',[5,1],4);
        Vec_DFF = imfilter(Vec_DFF,gausfilt);
        
        
        figure
        hold on
        for nn= 1:neurons
            shadedErrorBar([],squeeze(nanmean(Vec_DFF(:,:,nn),2)),...
                squeeze([ nanstd(Vec_DFF(:,:,nn),[],2)/sqrt(trials) * 1.96 ]));
        end
        Path_parts =  strsep(file_path,filesep);
        exp_path =  strrep(allfiles{file},'-','.');
        exp_parts  =  strsep(exp_path,filesep);
        title([Path_parts{end} ' '  exp_parts{end-1:end}])
        
        % smoothing
        fluoAllCorr = Output.fluoAllCorr;
        B_Vec = repmat(nanmean(fluoAllCorr(1:100,:,:)),[length(fluoAllCorr),1]);
        Vec_DFF2 = (fluoAllCorr -B_Vec)./B_Vec * 100;
        gausfilt = fspecial('gaussian',[11,1],4);
        Vec_DFF2 = imfilter(Vec_DFF2,gausfilt);
        
        FrameIdx = TimingInfo.FrameIdx;
        
        figure
        plot(Vec_DFF2)
        ax = axis;
        hold on
        for trial =  1:length(FrameIdx)
            plot([FrameIdx(trial,1),FrameIdx(trial,1)], [ax(3) ax(4)],'-r')
            plot([FrameIdx(trial,2),FrameIdx(trial,2)], [ax(3) ax(4)],'-r')
        end
        title([Path_parts{end} ' '  exp_parts{end-1:end}])
    end
    
   pause 
    
    
end
