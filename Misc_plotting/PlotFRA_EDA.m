function plotFRA_EDA(Out)
% FLO and DFF presorted by tone and levels


%% init
FreqLevels = unique(Out.FreqLevelOrder);
FreqLevels = sortrows(FreqLevels,{'Freqs','Levels'},{'Ascend','Descend'});
Freqs  = Out.FreqLevelOrder{:,1};
frqs   = length(Freqs);
uFreqs = unique(Freqs);
Levels = Out.FreqLevelOrder{:,2};
lvls   = length(Levels);
uLevels= unique(Levels);
H = height(FreqLevels);
DFF  = Out.DFF;
frames  = size(DFF,1);
trials  = size(DFF,2);
FLO = Out.FreqLevelOrder;
handles = Out.handles{1};
% parse a single psignalfile to get trial information

 % define behaviorally relevant timepoints
soundon = handles.PreStimSilence*handles.pfs;
soundoff = soundon + handles.PrimaryDuration * handles.pfs;
%% average response over trials

DFF2 = squeeze(nanmean(Out.DFF,2));        

%% sort by average  max response 

[~, DFF_order] = sort(max(DFF2));
DFF2 = DFF2(:,DFF_order);
DFF2 = DFF2(:,:);
%% plot average responses 
imagesc(DFF2')
FramesToSeconds
caxis([0 20])
%% Sort by time to peak and normalize responses to the max response 
[DFF2_max, DFF2_timing]  = max(DFF2);
% normalize
DFF3 = DFF2./DFF2_max; 


% sort
[DFF3_timing,DFF3_order]=  sort(DFF2_timing);
DFF3 = DFF3(:,DFF3_order);

% plot
figure
imagesc(DFF3')
caxis([0 1 ])
FramesToSeconds
title('Time To Max Response')


%% generate indices for different groups of data 
Noise_idx = FLO{:,2}== 0; % Noise at the same level as the pure tone
Tone_idx = FLO{:,2}== 99; % Pure tones
NoiseFRA_idx = ~Tone_idx; % All Noise 

% the same quantities but for the average matrix
mu_Noise_idx =FreqLevels{:,2} == 0;
mu_Tone_idx = FreqLevels{:,2}== 99;
mu_NoiseFRA_idx = ~mu_Tone_idx;

% Tones = DFF(:,Tone_idx,:);
% Noise = DFF(:,Noise_idx,:);
% NoiseFRA = DFF(:,NoiseFRA_idx,:);
% 

%% convert raw traces to mean response levels 
% converts into DFF_mu (Neuron X  unique Freq/levels X Frames) 

%remove all neruons with less than 2 percent DFF and all neurons with 
% a response of over 1000 on any single trial (noisy trials)
% (found empirically by reviewing DFF2 responses)
FluoroAll  = DFF(:,:,DFF_order);

% remove any trials with artifically high trials 
FluoroAll(repmat(max(FluoroAll)>700,150,1)) = nan; % trial clipping 
% restrict selection to neurons with average DFF above 2% and below 1000%
% - find better empirical upper bound -KS
DFFmin_idx = max(nanmean(FluoroAll,2 ))> 2; 
DFFmax_idx = max(max(FluoroAll)) < 500;
FluoroAll_idx = logical(squeeze(DFFmin_idx .* DFFmax_idx));
FluoroAll =   FluoroAll(:,:,FluoroAll_idx);

% preinitialize the loop
Neurons = size(FluoroAll,3);
DFF_mu = zeros(Neurons,H,frames);
% main loop
 for jj = 1:size(FluoroAll,3)
      neuron_ID = (jj);
      Fluoro = FluoroAll(:,:,neuron_ID);
     
      DFFtemp = zeros(height(FreqLevels),size(Fluoro,1));
      DFF_idx = 1 ;
      for kk = 1:height(FreqLevels)
          
              disp(' ')
             disp(strcat('Filtering Neuron: ', num2str(neuron_ID) ,';', ...
                  string(FreqLevels{kk,2}), ' dB; ',...
                  string(FreqLevels{kk,1}), ' Hz'))
              % find members of each Freq/Level to average
              idx =ismember(Freqs, FreqLevels{kk,1}) &...
                   ismember(Levels,FreqLevels{kk,2});
              
              n = sum(idx);
              if n > 0
%                   Ftemp = Fluoro(:,idx);
%                   F =  nanmean(Ftemp,2)';
%                   
%                   B =  nanmean(F(1:handles.PreStimSilence * handles.pfs));                    
%                   DF =( F - B)./B  ;
                 % DF = imfilter(DF, gausfilt);
                  DFFtemp(DFF_idx,:) = squeeze(...
                                      nanmean(DFF(:,idx,neuron_ID),2));
                    
                            
                
                  DFF_idx = DFF_idx + 1 ; 
              else
                  continue
              end
             
              
      end
     
      DFF_mu(jj,:,:) = DFFtemp;
        
     
 end

[df_by_level, df_by_level_timing ] = max(DFF_mu(:,:,:),[],3);

df_by_level_timing = df_by_level_timing + soundon - 1; % shift to correct idx

%% DFF cliping
% manually I've found that most trials with DFF higher than about 2-300 DFF
% are artifacts, consider using more robust methods to remove outliers -KS
%df_by_level(df_by_level>=500) = nan;

%% bf shift and Psignal Noise correlation changes 
% DFNoise = df_by_level(:,mu_Tone_idx);
% DFTones = df_by_level(:,mu_Noise_idx);
% DFNoiseFRA = df_by_level(:,mu_NoiseFRA_idx);
% 
% [BFNoise BFNoiseIdx] = max(DFNoise,[],2);
% [BFTone BFToneIdx] = max(DFTones,[],2);
% 
% % is there a higher response to noise or tones
% 
% BFTones_Noise = BFTone-BFNoise; 
% BFTones_Noise(BFTones_Noise > 1000 |BFTones_Noise < -1000 ) = nan;
% figure
% histogram(BFTones_Noise)
% nansum(BFTones_Noise>=0)
% 
% BFshift = BFToneIdx - BFNoiseIdx;
% figure
% histogram(BFshift,'Normalization','probability')
% unshifted = sum(abs(BFshift) <=2);
% % same analysis by quintile 
% q_idx = floor(linspace(0,Neurons,6)); % q_idx = quintile index
% for ii = 1:length(q_idx)-1
%   unshifted_quint(ii) =  sum(abs(BFshift(q_idx(ii)+1:q_idx(ii+1) )) <= 1);
%   fprintf(' %d %d \n' ,q_idx(ii)+1,q_idx(ii+1))
% end 
%sum(BF_Tones_Noise<1 & BF_Tones_Noise<-1)
% noise FRA 
    
% 
   
%% timing based classification
% look at noise responsive cells 


% %% Bandwith by Level
% bandwith_by_level = reshape(df_by_level,[Neurons,4,8] );
% 
% for ii =  1:Neurons
%     for jj = 1:4
%         temp_bandwith = df_by_level(ii,jj,:);
%         band_max = max(temp_bandwidth);
% %         
% %     end 
% % end 
% 
% 
% 
% % %% Significantly responding cells - This is likely depreciated 
% %                                     to update reshape to (150,10,4,8)
% %                                     frames X repeat X level X tone
% % Fluoro_Freqs = reshape(FluoroAll,150,8,40,[]);
% % Fluoro_Freqs = squeeze(mean(Fluoro_Freqs(60:100,:,:,:)));
% % Fluoro_Freqs = permute(Fluoro_Freqs,[2,1,3]);
% % sig_NT = zeros(Neurons,3);
% % sig_FRA = zeros(Neurons,3);
% % sig_CF  = zeros(Neurons,1);
% % Fluoro_Freqs_NT = Fluoro_Freqs(1:20,:,:); 
% % Fluoro_Freqs_FRA = Fluoro_Freqs(11:30,:,:); 
% % % FF_anova2 Freq response X Level Response X interaction
% % for ii = 1:Neurons
% %     sig_NT(ii,:)  = anova2(Fluoro_Freqs_NT(:,:,ii),10,'off'); 
% %     sig_FRA(ii,:) = anova2(Fluoro_Freqs_FRA(:,:,ii),10,'off');
% %     sig_CF_inf(ii,:) = anova1(Fluoro_Freqs(1:10,:,ii),[],'off');
% %     sig_CF_25(ii,:) = anova1(Fluoro_Freqs(11:21,:,ii),[],'off');
% %     sig_CF_15(ii,:) = anova1(Fluoro_Freqs(21:31,:,ii),[],'off');
% %     sig_CF_5(ii,:)  = anova1(Fluoro_Freqs(31:end,:,ii),[],'off');
% % end 
% % % turn into boolean at p =.05 significance Level
% %   NTsig_idx = sig_NT <.05;    
% %   FRAsig_idx = sig_FRA <.05;    
% % 
% %    alpha = .05
% %   for ii = 1:Neurons 
% %        if sig_CF_inf(ii)< alpha
% %            figure
% %           myFRA(uFreqs,uLevels,df_by_level(ii,:),ii)
% %      
% %        end
% %       if sig_CF_25(ii)< alpha
% %            figure
% %           myFRA(uFreqs,uLevels,df_by_level(ii,:),ii)
% %      
% %        end 
% %   end 
% %   
% %   
%   
%   for ii = 1:Neurons 
%       if  NTsig_idx(ii,1)
%           figure
%           myFRA(uFreqs,uLevels,df_by_level(ii,:),ii)
%           disp(' Freq selective neuron')
%           pause
%       end 
%       
%       
%       if NTsig_idx(ii,3)
%           figure
%           myFRA(uFreqs,uLevels,df_by_level(ii,:),ii)
%           disp('Freq&Sound level selective')
%           pause
%       end 
% 
%   end 
% 
%   
%   NT_level_idx = sig_NT <.05;
%   NT_level_idx = NT_level_idx(:,2);
%   
%   
%   
%   
%   FF_levels = Fluoro_Freqs_NT(:,:,NTsig_idx(:,2));
%   
% FF_Tone_larger_idx = squeeze( sum(sum((FF_levels(1:10,:,:))))>  sum(sum((FF_levels(2:20,:,:)))) )
 %%  Psignal Correlations 
 % Assuming All Neurons were V-Shaped or tuned in a similar way. one would
 % expect that signal correlation would be highest in the INF SNR condition
 % and decrease as noise was added This code tests that assumption
 
 % init
%  s_cov = zeros(Neurons,Neurons);
%  s_corr= zeros(Neurons,Neurons,4);
%  idx = [1:10:41];
%  DFF_index = max(max(Fluoro_Freqs))>100;               % index %true(Neurons,1)== All neurons
%  reshape_num =sum(DFF_index)      % # Neurons; 
%  % main loop
%  for ii = 1:4
%      FFtemp = Fluoro_Freqs(idx(ii):idx(ii+1)-1,:,DFF_index);
%      FFtemp = reshape(FFtemp,[], reshape_num ) ;
%      s_cov = cov(FFtemp);
%      s_temp = getCorrFromCov(s_cov);
%      s_corr_100(:,:,ii) = s_temp;
%      s_mean(ii)= nanmean(s_temp(:));
%      s_std(ii) = nanstd(s_temp(:));
%  end 
%   clear FFtemp
%   clear s_temp 
%   figure
%   hold on 
%   errorbar(1:4,s_mean+.05,s_std/2,'.')
%   bar(1:4,s_mean+.05)
%    
%   figure
%   for ii = 1:4
%       subplot(1,4,ii) 
%       corr_temp = s_corr_100(:,:,ii);
%       histogram(corr_temp,20,'Normalization','probability')
%   end
%   
%   clear corr_temp 
  
 %%  optional debugging/verification
 % this plot should show that taking the grand average either by all trials
 % or an average of each average response both result in the same plot
%  subplot(1,2,1)
%  plot(squeeze(mean(DFF_mu(1,:,:))))
%  subplot(1,2,2)
%  plot(squeeze(mean(DFF(:,:,1)')))



% create average sound evoked response by level
 df_by_level =  max(DFF_mu(:,:,soundon:soundoff),[],3) ;

 
 % create grand average timecourse and filter it 
 DFF_grand_avg = squeeze(nanmean(DFF_mu,2));
 
 gausfilt = fspecial('gaussian',[20,1],5);
 for ii = 1:Neurons 
 DFF_grand_avg(ii,:) = imfilter(DFF_grand_avg(ii,:),gausfilt);
 end 

for ii = 1:Neurons
      
    % extract a single neurons Responses
   
   
   
  % plot the responses 
    figure
    set(gcf,'Position',[5   129   648   837])
    subplot(1,2,1)
    mu = DFF_grand_avg(ii,:);
    hold on 
    plot(mu,'b')
    %plot(nanmean(Vec_DFF(:,:,nn),2),'k','LineWidth',2)
    aa = axis;
    plot([aa(1), aa(2)], [0 0 ] ,'--k')
    plot([30 30], [aa(3),aa(4)],'--r')
    plot([60 60], [aa(3),aa(4)],'--g')
    plot([90 90], [aa(3),aa(4)],'--g')
    plot([120 120], [aa(3),aa(4)],'--r')
    

    subplot(1,2,2)
    myFRA(uFreqs,uLevels,df_by_level(ii,:),ii);
    % yticklabels({'inf','20' ,'+10','00'})
    
    

end
end 




function FramesToSeconds 
ax = gcf;
figure(ax) 
xticks([30,60,90,120])
xticklabels({'1','2','3','4'})
xlabel('Time (sec)')
ylabel('Neurons')
end 
