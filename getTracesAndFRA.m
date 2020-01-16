function[df_by_level,DFF2,DFF3] = getTracesAndFRA(Raw)
% this function is a scripted function for Patrick to play with FRAs and
% Average responses of Cells -KS 
%% Input  
% make the folder that is function is contained in your working directory
% as well as have the associated psignal file 


%% Output variables 
% DFF2 is a matrix of each neurons average response across Levels 
% df_by_level is matrix of each neurons maximum response for each 
%tone/level combination 

%% init
FreqLevels = unique(Raw.FreqLevelOrder);
FreqLevels = sortrows(FreqLevels,{'Freqs','Levels'},{'Ascend','Descend'});
Freqs  = Raw.FreqLevelOrder{:,1};
frqs   = length(Freqs);
uFreqs = unique(Freqs);
Levels = Raw.FreqLevelOrder{:,2};
lvls   = length(Levels);
uLevels= unique(Levels);
H = height(FreqLevels);
DFF  = Raw.DFF;
frames  = size(DFF,1);
trials  = size(DFF,2);
FLO = Raw.FreqLevelOrder;

% parse a single psignalfile to get trial information
if iscell(Raw.DataDirs)
    
    Main_path = Raw.DataDirs{1};
else 
    Main_path = Raw.DataDirs
end 

dir_t = dir([Main_path] );
dir_t = {dir_t.name};
f_bl= ~ cellfun(@isempty,(regexp( dir_t  ,'_Phys_')));

Psignal_file= dir_t{f_bl};

% Psignal Handling

handles = WF_getPsignalInfo(fullfile(Main_path, Psignal_file));

 % define behaviorally relevant timepoints
soundon = handles.PreStimSilence*handles.pfs;
soundoff = soundon + handles.PrimaryDuration * handles.pfs;

%% average response over trials

DFF2 = squeeze(mean(Raw.DFF,2));        


%% sort by average  max response 

[~, DFF_order] = sort(max(DFF2));
DFF2 = DFF2(:,DFF_order);
DFF2 = DFF2;





%% plot average responses 
figure
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

% s =elect only active neurons 
active_idx = Raw.Active{:,2} > 0 ;
active_idx = active_idx(DFF3_order);
DFF3 = DFF3(:,active_idx);

% plot
figure
imagesc(DFF3')
caxis([0 1 ])
FramesToSeconds
title('Time To Max Response')



%% convert raw traces to mean response levels 
% converts into DFF_mu (Neuron X  unique Freq/levels X Frames) 

%remove all neruons with less than 2 percent DFF and all neurons with 
% a response of over 1000 on any single trial (noisy trials)
% (found empirically by reviewing DFF2 responses)
FluoroAll  = DFF(:,:,DFF_order);

% remove any trials with artifically high trials 
% FluoroAll(repmat(max(FluoroAll)>700,150,1)) = nan; % trial clipping 
% % restrict selection to neurons with average DFF above 2% and below 1000%
% % - find better empirical upper bound -KS
% DFFmin_idx = max(mean(FluoroAll,2 ))> 2; 
% DFFmax_idx = max(max(FluoroAll)) < 500;
% FluoroAll_idx = logical(squeeze(DFFmin_idx .* DFFmax_idx));
% FluoroAll =   FluoroAll(:,:,FluoroAll_idx);
% 
% % restrict corresponding DFF2 traces
% DFF2 = DFF2(FluoroAll_idx,:);

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

                  DFFtemp(DFF_idx,:) = squeeze(...
                                      nanmedian(FluoroAll(:,idx,neuron_ID),2));
                    
                            
                
                  DFF_idx = DFF_idx + 1 ; 
              else
                  continue
              end
             
              
      end
     
      DFF_mu(jj,:,:) = DFFtemp;
        
     
 end

[df_by_level, ~ ] = max(DFF_mu(:,:,:),[],3);
%[df_by_level_tones] = max(DFF_mu(:,:,60:120),[],3);


