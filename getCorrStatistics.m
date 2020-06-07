function [corr_mu,corr_std,corr_CI,corr_sig] =...
                 getCorrStatistics(corrs,name)
    % this function takes the correlations and flattens them and extracts
    % the relevant statistics 
    if iscell(corrs)
        corr_flat = cell2mat(cellfun(@(x)x(:),corrs,'UniformOutput',0));
    else 
        corr_flat = corrs(:);
    end 
corr_mu = nanmean(corr_flat);
corr_std = nanstd(corr_flat); 
corr_CI = corr_std ./ sqrt(length(corr_flat))  * 1.96;  % 95 percent Confidence interval

[p ~,stats] = anova1(corr_flat,[],'off');
try
corr_sig = multcompare(stats,[],'off');
catch 
    corr_sig =[];
end 