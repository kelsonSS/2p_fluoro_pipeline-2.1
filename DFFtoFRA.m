function df_by_level,Sig_by_level = DFFtoFRA(In,Sig_flag);
% this function takes a data container from Fluoro_to_table and takes the DFF and Trial
% information therein and outputs the average response to each frequency
% level combination

Sig_by_level = 0;

DFF = In.DFF;

FreqLevels = unique(In.FreqLevelOrder);
Freqs =In.FreqLevelOrder{:,1};
Levels= In.FreqLevelOrder{:,2};
uLevels= unique(Levels);
uFreqs = unique(Freqs);
try
active_list = In.Active.Activity;
catch
 active_list = In.active.Activity;   
end 

handles = In.handles;
if length(handles) >1
    handles = handles(1);
end 
soundon = handles.PreStimSilence * handles.pfs;
soundoff = soundon + handles.PrimaryDuration + handles.pfs; 

for jj = 1:size(DFF,3)
      neuron_ID = (jj);
      Fluoro = DFF(:,:,neuron_ID);
     
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
[df_by_level_baseline] = max(DFF_mu(1:30,:,:),[],3);

    








end 
 

