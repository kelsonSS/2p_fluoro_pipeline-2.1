function Out = getFluoroHitsMiss(TNBehavior,out_folder)
    
if ~exist('out_folder','var')
    to_plot = 0
else
    to_plot = 1
end 

noise_level = 50;
detectLevelSNR = [-10, 0 ,10,20];

Out = CreateOutputStructure(detectLevelSNR);


detectLevels = detectLevelSNR  + noise_level;    

    for expt_idx =1:length(TNBehavior)
        
        
        expt= TNBehavior{expt_idx};
        handles = expt.handles{1};
        uLevels = handles.uLevels - 50;
        
        F = expt.DFF_Z(:,:,expt.Clean_idx & expt.active{:,2}>0);
        
        hits_idx = handles.Hits';
        miss_idx = handles.Miss';
        early_idx = handles.Early';
        
        all_lvl_idx = any(handles.Levels - noise_level == [-10 0 10 20],2);
        
        all_hits_temp = squeeze(nanmean(F(:,all_lvl_idx & hits_idx,:),2));
        all_miss_temp = squeeze(nanmean(F(:,all_lvl_idx & miss_idx,:),2));
        all_early_temp = squeeze(nanmean(F(:,all_lvl_idx & early_idx,:),2));
        
        
        
        Out.All.Hits = cat(2,Out.All.Hits,all_hits_temp);
        Out.All.Miss = cat(2,Out.All.Miss,all_miss_temp);
        Out.All.Early = cat(2,Out.All.Early,all_early_temp);

        
        for lvl = 1:length(uLevels)
            uLevel = uLevels(lvl);
            if any(uLevel == detectLevelSNR)
                
                lvl_idx = (handles.Levels - noise_level) == uLevel;
                hits_temp = squeeze(nanmean(F(:,lvl_idx & hits_idx,:),2));
                miss_temp = squeeze(nanmean(F(:,lvl_idx & miss_idx,:),2));
                early_temp = squeeze(nanmean(F(:,lvl_idx & early_idx,:),2));
             
                
                
                
                snr_name = sprintf('SNR_%d',uLevel);
                snr_name = strrep(snr_name,'-','minus_');
               
                Out.(snr_name).Hits =cat(2,Out.(snr_name).Hits,hits_temp);
                Out.(snr_name).Miss =cat(2,Out.(snr_name).Miss,miss_temp);
                Out.(snr_name).Early = cat(2,Out.(snr_name).Early,early_temp);

            end 
        end 
    end 
    
    if to_plot
        PlotResults(Out,out_folder)
    end 

    
end 



function Out =CreateOutputStructure(levels)

levels = [num2cell(levels),'All'];

for ii = 1:length(levels)
if isnumeric(levels{ii})
    snr_name = sprintf('SNR_%d',levels{ii});
    snr_name = strrep(snr_name,'-','minus_');
else 
    snr_name = levels{ii};
end 

Out.(snr_name).Hits = [];
Out.(snr_name).Miss = [];
Out.(snr_name).Early = [];
end 


end 

function PlotResults(Out,out_path)

f = fieldnames(Out);

for lvl =1:length(f)
    
figure; hold on
shadedErrorBar([],nanmean(Out.(f{lvl}).Hits,2),nanstd(Out.(f{lvl}).Hits,[],2)/ sqrt(length(Out.(f{lvl}).Hits)) * 1.96, 'b')
shadedErrorBar([],nanmean(Out.(f{lvl}).Miss,2),nanstd(Out.(f{lvl}).Miss,[],2)/ sqrt(length(Out.(f{lvl}).Miss)) * 1.96 )
shadedErrorBar([],nanmean(Out.(f{lvl}).Early,2),nanstd(Out.(f{lvl}).Early,[],2)/ sqrt(length(Out.(f{lvl}).Early)) * 1.96,'r')
title_name = sprintf('%s',f{lvl});
title(title_name,'Interpreter','none');
print(fullfile(out_path,[title_name '.pdf']),'-dpdf')

figure
shadedErrorBar([],nanmean(Out.(f{lvl}).Hits-Out.(f{lvl}).Miss,2),nanstd(Out.(f{lvl}).Hits-Out.(f{lvl}).Miss,[],2)/ sqrt(length(Out.(f{lvl}).Hits)) * 1.96)
title_name =sprintf('%s: Hit-miss',f{lvl});
title(title_name ,'Interpreter','none');


print(fullfile(out_path, sprintf('%s_HitMissDiff.pdf',f{lvl}) ),'-dpdf' ) 



end 
end 
