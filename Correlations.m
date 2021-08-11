function Out = Correlations(TN,normalized)

if ~exist('normalized','var')
    normalized = 0
end 

%Cells = TN.CellID; 
if normalized 
    DFF = TN.DFF_Z;
else 
DFF   = TN.DFF;
end 
active = TN.active.Activity;
FreqLevelOrder = TN.FreqLevelOrder;

L_num = length(unique(FreqLevelOrder.Levels));
F_num = length(unique(FreqLevelOrder.Freqs));
e_ls = TN.experiment_list;
rep_mode = 5; % trials per unique condition
n_frames = size(DFF,1);
n_trials =size(DFF,2);
handles = TN.handles{1};

soundon = handles.PreStimSilence * 30; % everything is 30FPS 
soundoff = soundon + handles.PrimaryDuration * 30; 

% if there is noise
if handles.BackgroundNoise(1) ~= -99
    soundon = min(soundon,handles.BackgroundNoise(2) * 30);
    soundoff = max(soundoff,handles.BackgroundNoise(3) * 30);
end 

 % rehape Frame X Repeat X Level X Frequency X Neuron 


Expt_num = length(TN.DataDirs);
for Expt = 1:Expt_num
    % unpack cells for this experiment and filter to only responsive cells      
   %  cells = TN.CellID{Expt}; could be used for spatial locations of each
   %  cell
   if exist('L','var')
       fprintf(repmat('\b',1,progress_text_length)) % clear previous line
   end 
   
   % print progress 
   progress_text= sprintf('analyzing expt %d of %d',Expt, Expt_num);
    progress_text_length = length(progress_text);
    fprintf(progress_text);
   
    DFF_t = DFF(:,:,e_ls == Expt);
    a_t  = (active(e_ls == Expt) >= 1); % <= P.05 significance value 
   DFF_t = DFF_t(:,:,a_t);
    n_neurons = size(DFF_t,3);
%   if isfield( TN,'handles')
%         FLO{:,1} = TN.handles{Expt}.Freqs;
%         FLO{:,2} = TN.handles{Expt}.Levels;
%         uFreqs = unique(FLO{:,1});
%         uLevels = unique(FLO{:,2});
%         L = length(uLevels);
        
%   else 
       FLO = FreqLevelOrder;
        Levels = sort(unique(FLO{:,2}),'descend');
        L = length(uLevels);
        uFreqs = unique(FreqLevelOrder.Freqs);
%   end 
  
    % continue if  size  = 1 
    if n_neurons <= 1 
         Corr{Expt} = nan;
    for Lvl = 1:L
    LCorr{Expt,Lvl} = nan;
    NCorr{Expt,Lvl}  = nan;
    end 
    continue
    end 
        
   
    
    
     % rehape Frame X Repeat X Level X Frequency X Neuron 
      for freq = 1:F_num
        for lvl = 1:L_num
         
             fprintf('Level %d Freq %d \n', uLevels(lvl),uFreqs(freq))
           FL_idx=  FLO{:,1} == uFreqs(freq) &... 
           FLO{:,2} == uLevels(lvl);
            
            Vec_DFF_Temp = DFF_t(:,FL_idx,:);
            % in case trials have different reps
            try
            Vec_DFF_Temp = Vec_DFF_Temp(:,1:rep_mode,:);
            catch
            end
          %  plot(Vec_DFF_Temp(:,:,1))
            % baseline = squeeze(nanmean(Vec_DFF_Temp(1:30,:,:),2));
            % get mean response to this freq/Level to tones 
            temp_mu =  squeeze(nanmean(nanmean(Vec_DFF_Temp(soundon:soundoff,:,:),2)));
            DFF_mu(lvl,freq,:) =temp_mu;
            DFF_n_temp = squeeze(...
                       nanmean( Vec_DFF_Temp(soundon:soundoff,:,:),1));
                    
            DFF_n_temp = DFF_n_temp - squeeze(...
                                               repmat( DFF_mu(lvl,freq,:),...
                                                    [rep_mode,1] ) );      
            DFF_n2(lvl,freq,:,:) = DFF_n_temp;
            DFF_c (:,:,freq,lvl,:) = Vec_DFF_Temp;
              
              clear Vec_DFF_Temp
      
        end
    end 
     
     
   % DFF_c =  reshape(DFF_t, [frames,[],L_num,F_num,N_t]);
    %DFF_mu = squeeze(nanmean(nanmean(DFF_c(60:120,:,:,:,:),2))); % average over stim frames and repeats
    % change to column major form 
    DF_flat = reshape(DFF_mu,[L_num *F_num ,size(DFF_mu,3)]); 
    
 
    
    
    % subtract the mean for each sound/level to find residuals
   % DFF_n = DFF_n2 - permute(repmat(DFF_mu,1,1,1,10),...
   %                                  [4 2 1 3 ]);
                                 
   % DFF_n2 = reshape(DFF_n2, 10 * 7 , 4, [] );           
    
  
    
    
    
    
    Corr{Expt} = getCorrFromCov(cov(DF_flat));
    clear DF_flat
    for Lvl = 1:size(DFF_mu,1);
    LCorr{Expt,Lvl} = getCorrFromCov(cov(squeeze(DFF_mu(Lvl,:,:))));
    %NCorr{Expt,Lvl}  = getCorrFromCov(cov(squeeze(DFF_n(:,Lvl,:))));
    
    NCorr2{Expt,Lvl} = getCorrFromCov( cov(...
                reshape( DFF_n2(Lvl,:,:,:),...
                         [] , n_neurons) ) );
    end 
 
 
      clear DFF_mu
      clear DFF_n2
      clear DFF_c
      
end 
%% Correlation analysis
%analysis

disp('Analyzing Noise Statistics')
%[N_mu,N_std,N_CI,N_sig] = getCorrStatistics(NCorr,'Noise');
[N_mu,N_std,N_CI,N_sig,N_corr_flat] = getCorrStatistics(NCorr2,'Noise2'); 
disp('Analyzing Signal Statistics')
[L_mu,L_std,L_CI,L_sig,L_corr_flat] = getCorrStatistics(LCorr,'Signal'); 
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
% figure;cdfplot(L_corr_flat);title('Signal Correlation')
% figure;cdfplot(N_corr_flat);title('Noise Correlation')

% plotting- Ncorr Bars
Bars_by_level(L_mu,L_CI,'Signal')
Bars_Overall(L_mu,L_CI,'Signal')
Bars_by_level(N_mu,N_CI,'Noise')
Bars_Overall(N_mu,N_CI,'Noise')
% Bars_by_level(N_mu2,N_CI2,'Noise2')
%% packaging
    Out.Corr = Corr;
    Out.LCorr = LCorr;
    Out.NCorr = NCorr2;
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


function [corr_mu,corr_std,corr_CI,corr_sig,corr_flat] =...
                 getCorrStatistics(corrs,name)
    % this function takes the correlations and flattens them and extracts
    % the relevant statistics 
corr_flat = cell2mat(cellfun(@(x)x(:),corrs,'UniformOutput',0));

corr_mu = nanmean(corr_flat);
corr_std = nanstd(corr_flat); 
corr_CI = corr_std ./ sqrt(length(corr_flat))  * 1.96;  % 95 percent Confidence interval

[p ~,stats] = anova1(corr_flat,[],'off');
if size(corr_flat,2) >1
    corr_sig = multcompare(stats,[],'off')
else 
    corr_sig = p;
end 
% figure
% hold on 
% for ii = 1:4
% cdfplot(corr_flat(:,ii))
% end
% title(sprintf('%s Correlation by Level',name))
% 



