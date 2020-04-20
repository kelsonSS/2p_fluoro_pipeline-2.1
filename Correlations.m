function Out = Correlations(TN)

%Cells = TN.CellID; 
DFF   = TN.DFF;
active = TN.active.Activity;
FLO = TN.FreqLevelOrder;
L_num = unique(FLO.Levels);
F_num = unique(FLO.Freqs);
e_ls = TN.experiment_list;

frames = size(DFF,1);
trials =size(DFF,2);
neurons = size(DFF,3);
 % rehape Frame X Repeat X Level X Frequency X Neuron 



for Expt = 1:length(TN.DataDirs)
    % unpack cells for this experiment and filter to only responsive cells      
   %  cells = TN.CellID{Expt}; could be used for spatial locations of each
   %  cell
    DFF_t = DFF(:,:,e_ls == Expt);
    a_t  = (active(e_ls == Expt) > 2); % <= P.01 significance value 
   DFF_T = DFF_t(:,:,a_t);
    N_t = size(DFF_t,3);
   
    % continue if  size  = 1 
    if N_t <= 1 
         Corr{Expt} = nan;
    for Lvl = 1:4
    LCorr{Expt,Lvl} = nan;
    NCorr{Expt,Lvl}  = nan;
    end 
    continue
    end 
        
   
    
     % rehape Frame X Repeat X Level X Frequency X Neuron 
    DFF_c =  reshape(DFF_t, [frames,10,L_num,F_num,N_t]);
    DFF_mu = squeeze(nanmean(nanmean(DFF_c(60:120,:,:,:,:),2))); % average over stim frames and repeats
    % change to column major form 
    DF_flat = reshape(DFF_mu,[32,size(DFF_mu,3)]); 
    
    % arrange data for Noise Correlations 
    DFF_n = permute(DFF_c,[1,2,4,3,5]);
    DFF_n = reshape(DFF_n,[(150 * 10 * 8) , 4 , N_t]);
    
    % standard form of noise corr where we subtract residuals from tone
    % period 
    DFF_n2 = squeeze(nanmean(DFF_c(60:120,:,:,:,:)));
    DFF_n2 = permute(DFF_n2,[1,3,2,4]);
    
    
    
    % subtract the mean for each sound/level to find residuals
    DFF_n2 = DFF_n2 - permute(repmat(DFF_mu,1,1,1,10),...
                                     [4 2 1 3 ]);
                                 
    DFF_n2 = reshape(DFF_n2, 10 * 8 , 4, [] );           
    

    
    
    
    
    Corr{Expt} = getCorrFromCov(cov(DF_flat));
    for Lvl = 1:4
    LCorr{Expt,Lvl} = getCorrFromCov(cov(squeeze(DFF_mu(Lvl,:,:))));
    NCorr{Expt,Lvl}  = getCorrFromCov(cov(squeeze(DFF_n(:,Lvl,:))));
    NCorr2{Expt,Lvl} = getCorrFromCov(cov(squeeze(DFF_n2(:,Lvl,:))));
    end 
 
end 
      
%% Correlation analysis
%analysis
[N_mu,N_std,N_CI,N_sig] = getCorrStatistics(NCorr,'Signal');
[N_mu2,N_std2,N_CI2,N_sig2] = getCorrStatistics(NCorr2,'Noise'); 
[L_mu,L_std,L_CI,L_sig] = getCorrStatistics(LCorr,'Noise2'); 
%  Older version - Delete if above is running fine 
%  LCorr_flat = cell2mat(cellfun(@(x)x(:),LCorr,'UniformOutput',0));
%  NCorr_flat = cell2mat(cellfun(@(x)x(:),NCorr,'UniformOutput',0));
%  NCorr2_flat = cell2mat(cellfun(@(x)x(:),NCorr2,'UniformOutput',0));
%  
% N_mu = nanmean(NCorr_flat);
% N_std = nanstd(NCorr_flat); 
% CI = N_std ./ sqrt(length(NCorr_flat))  * 1.96;  % 95 percent Confidence interval
% 
% N2_mu = nanmean(NCorr2_flat);
% N2_std = nanstd(NCorr2_flat); 
% CI2 = N2_std ./ sqrt(length(NCorr2_flat))  * 1.96;  % 95 percen
% 
% 
% [p ~,stats] = anova1(LCorr_flat,[],'off');
% Sig_S = multcompare(stats,[],'off'); 
% 
% [p ~,stats] = anova1(NCorr_flat,[],'off');
% Sig_N = multcompare(stats,[],'off'); 
% 
% [p ~,stats] = anova1(NCorr2_flat,[],'off');
% Sig_N2 = multcompare(stats,[],'off'); 


%% plotting - Signal and noise CDFs

%


% plotting- Ncorr Bars
Bars_by_level(L_mu,L_CI,'Signal')
Bars_by_level(N_mu,N_CI,'Noise')
Bars_by_level(N_mu2,N_CI2,'Noise2')
%% packaging
    Out.Corr = Corr;
    Out.LCorr = LCorr;
    Out.NCorr = NCorr;
    Out.Signal_Sig = L_sig;
    Out.Noise_Sig  = N_sig;
    Out.Signal_Stats.mean = L_mu;
    Out.Signal_Stats.std = L_std;
    Out.Noise_Stats.mean = N_mu;
    Out.Noise_Stats.std = N_std;
    

    
    
function Bars_by_level(mu,CI,name)
       % takes a mean and 95% confidence interval and plots them 
        % plotting- Ncorr Bars
figure
errorbar(mu,CI, '.')
hold on
bar(mu,'k')
hold off
title(sprintf('%s Correlations by level',name))
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',{'INF','+20','+10','0'})


function [corr_mu,corr_std,corr_CI,corr_sig] =...
                 getCorrStatistics(corrs,name)
    % this function takes the correlations and flattens them and extracts
    % the relevant statistics 
corr_flat = cell2mat(cellfun(@(x)x(:),corrs,'UniformOutput',0));

corr_mu = nanmean(corr_flat);
corr_std = nanstd(corr_flat); 
corr_CI = corr_std ./ sqrt(length(corr_flat))  * 1.96;  % 95 percent Confidence interval

[p ~,stats] = anova1(corr_flat,[],'off');
corr_sig = multcompare(stats,[],'off')

% figure
% hold on 
% for ii = 1:4
% cdfplot(corr_flat(:,ii))
% end
% title(sprintf('%s Correlation by Level',name))
% 



            