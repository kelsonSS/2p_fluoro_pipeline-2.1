function Out = CorrelationsByAnimalBehavior(TNBehavior,uLevels,Timing)

normalized = 0; % zero for DFF 1 for DFF_Z

findLevels_flg = false;
if ~exist('uLevels','var') || isempty(uLevels)
    uLevels = [];
    findLevels_flg = true;
end 


if ~exist('Timing','var')
    Timing = 'Tone';
end 



Expt_num = length(TNBehavior);
Expt_ID = 1;
for Expt = 1:Expt_num
    TN = TNBehavior{Expt};
    if normalized
        DFF = TN.DFF_Z;
    else
        DFF   = TN.DFF;
    end
    active = TN.active.Activity;
    
        FLO = TN.FreqLevelOrder;
        
        if findLevels_flg
            uLevels = sort(unique(FLO.Levels),'descend');
        end 
        L = length(uLevels);
        uFreqs = unique(FLO.Freqs);
        uFreqs = uFreqs(uFreqs < 60000);
    
    rep_mode = 10; % trials per unique condition
    
    lvls_to_check = any(FLO{:,2} == uLevels,2);
    counts = histcounts2(FLO{lvls_to_check,1},FLO{lvls_to_check,2});
    counts(counts == 0 ) =[];
    rep_mode = min(counts);
    
    
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


    % if specifically looking at prestim activity
    if strcmp(Timing,'PreStim')
        soundon = 1;
        soundoff = soundon + handles.PreStimSilence * 30;
    end 
    
    % unpack cells for this experiment and filter to only responsive cells      
   %  cells = TN.CellID{Expt}; could be used for spatial locations of each
   %  cell
   if exist('progress_text_length','var')
       fprintf(repmat('\b',1,progress_text_length)) % clear previous line
   end 
   
   % print progress 
   progress_text= sprintf('analyzing expt %d of %d',Expt, Expt_num);
    progress_text_length = length(progress_text);
    fprintf(progress_text);
   
    a_t  = active >= 1; % <= P.05 significance value 
   DFF_t = DFF(:,:,a_t);
    n_neurons = size(DFF_t,3);

   
% 
    % continue if  size  = 20 
     if any([n_trials,n_neurons] <= 1) 
          Corr{Expt,1} = nan;
          NCorrAll{Expt,1} = nan;
          
       continue
     end 
        
   
    
      F_num = length(uFreqs);
      
      L_num = length(uLevels);
      
     % rehape Frame X Repeat X Level X Frequency X Neuron 
      for freq = 1:F_num
        for lvl = 1:L_num
         
              %fprintf('Level %d Freq %d \n', uLevels(lvl),uFreqs(freq))
           FL_idx=  FLO{:,1} == uFreqs(freq) &... 
           FLO{:,2} == uLevels(lvl);
            
            Vec_DFF_Temp = DFF_t(:,FL_idx,:);
            % in case trials have different reps
            try
            Vec_DFF_Temp = Vec_DFF_Temp(:,1:rep_mode,:);
            catch
                Vec_DFF_Temp = nan(size(Vec_DFF_Temp,1),rep_mode,size(Vec_DFF_Temp,3));
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
   % DFF_n = DFF_n2 - permute(repmat(DFF_mu,1,1,1,rep_mode),...
   %                                  [1 2 4 3 ]);
                                 

    
    
    
    
    Corr{Expt_ID,1} = getCorrFromCov(cov(DF_flat));
    
  
    
    NCorrAll{Expt_ID,1} = getCorrFromCov( cov(...
                 reshape( DFF_n2(:,:,:,:),...
                          [] , n_neurons) ) );
    
    clear DF_flat
     for Lvl = 1:size(DFF_mu,1)
 
   %     LCorr{Expt_ID,Lvl} = getCorrFromCov(cov(squeeze(DFF_mu(Lvl,:,:))));
        
        cov_mat = cov( reshape( DFF_n2(Lvl,:,:,:),[],n_neurons) );
    
        if all(isnan(cov_mat(:))) || length(cov_mat) == 1 % if we cannot find covarience return NAs
            NCorr2{Expt_ID,Lvl} = nan(n_neurons);
        else 
            
            NCorr2{Expt_ID,Lvl} = getCorrFromCov( cov_mat );
       
        end 
        
     end
     
     Expt_ID = Expt_ID + 1 ;
      clear DFF_mu
      clear DFF_n2
      clear DFF_c
      
end 
%% Correlation analysis
%analysis

disp('Analyzing Noise Statistics')
[~,~,~,~,~,NCorrAnimal] = getCorrStatistics(NCorr2,'Noise');
[~,~,~,~,~,NCorrAllAnimal] = getCorrStatistics(NCorrAll,'Noise2'); 
%disp('Analyzing Signal Statistics')
[~,~,~,~,~,CorrAnimal] = getCorrStatistics(Corr,'Signal'); 
%[~,~,~,~,~,LCorrAnimal] = getCorrStatistics(LCorr,'Signal');
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
%Bars_by_level(L_mu,L_CI,'Signal');AddScatter(L_mu_expt)
%Bars_Overall(L_mu,L_CI,'Signal')
%Bars_by_level(N_mu,N_CI,'Noise');AddScatter(N_mu_expt)
%Bars_Overall(N_mu,N_CI,'Noise')
% Bars_by_level(N_mu2,N_CI2,'Noise2')


%% packaging
    Out.Corr = Corr;
   % Out.LCorr = LCorr;
   % Out.LCorrAnimal = LCorrAnimal;
    Out.NCorrAll = NCorrAll;
    Out.NCorrAnimal = NCorrAllAnimal;
    Out.NCorr = NCorr2;
    Out.NCorrAnimal = NCorrAnimal;
    %Out.NCorr = NCorr2;
    %Out.NCorrAnimal = N_mu_expt;
    %%Out.Signal_Sig = L_sig;
    %Out.Noise_Sig  = N_sig;
    %Out.Signal_Stats.mean = L_mu;
    %Out.Signal_Stats.std = L_std;
    %Out.Noise_Stats.mean = N_mu;
    %Out.Noise_Stats.std = N_std;
    %Out.LCorrDiff = L_diff;
    %Out.LCorrDiffAnimal = L_diff_expt;
    %Out.NCorrDiff = N_diff;
    %Out.NCorrDiffAnimal = N_diff_expt;
    



        
        
    
function Bars_by_level(mu,CI,name)
       % takes a mean and 95% confidence interval and plots them 
        % plotting- Ncorr Bars
figure
bar(mu)
hold on
errorbar(mu,CI, 'k.')
hold off
title(sprintf('%s Correlations by level',name))
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',{'INF','+20','+10','0'})


function AddScatter(mu_expt)
        % adds scatterplot to the current figure
        % currently only used after Bars_by_level
        hold on 
        levels = size(mu_expt,2);
        n_expts = size(mu_expt,1);
        
        timing = repmat( [1:levels],n_expts,1);
        
        scatter(timing(:),mu_expt(:),'k.' )
        
        
        



function [corr_mu,corr_std,corr_CI,corr_sig,corr_flat,corr_mu_expt,...
          corr_diff,corr_diff_expt] =...
                 getCorrStatistics(corrs,name)
    % this function takes the correlations and flattens them and extracts
    % the relevant statistics 

CorrStats = struct()   ; 
corr_flat = cellfun(@(x)x(:),corrs,'UniformOutput',0);


corr_mu_expt =   cellfun(@nanmean ,corr_flat);
corr_mu = nanmean(corr_mu_expt);
corr_std =  nanstd(corr_mu_expt); 
corr_CI = corr_std ./ sqrt(length(corr_flat))  * 1.96;  % 95 percent Confidence interval

[p ~,stats] = anova1(cell2mat(corr_flat),[],'off');
if size(corr_flat,2) >1
    corr_sig = multcompare(stats,[],'off');
else 
    corr_sig = p;
end 

% see how correlation changes in comparision to tone level
corr_diff = {};
for ii = 1:size(corrs,1)
    corr_diff(ii,:) = cellfun( @(x) x - corrs{ii,1},corrs(ii,:),'UniformOutput',false);
end 

corr_diff_expt =  cellfun(@(x) nanmean(x(:)),corr_diff);


% figure
% hold on 
% for ii = 1:4
% cdfplot(corr_flat(:,ii))
% end
% title(sprintf('%s Correlation by Level',name))
% 



