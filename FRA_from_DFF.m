function df_by_level = FRA_from_DFF(FluoroAll, FreqLevels) 

soundon  = 60;
soundoff = 90;
Freqs = FreqLevels{:,1};
Levels = FreqLevels{:,2};
FreqLevels = sortrows(unique(FreqLevels),[1 -2]);

NNeurons = size(FluoroAll,3);
for ii = 1:NNeurons
    N = FluoroAll(:,:,ii)';
[~,p]=ttest2(mean(N(:,1:soundon),2),mean(N(:,soundon+1:soundoff),2));
  active(ii,1) = p;
end 
   % find significantly active neurons  number indicates the number of
   % stars that would be added to the figure. Zero indicates inactivity
   % while anything above 1 is active
   %
   active(active <.001) = 3;  
   active(active <.01)  = 2;  
   active(active <=.05) = 1;  
   active(active < 1  ) = 0;     


   neuron_ID = find(active >2);
   
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
                                      nanmedian(FluoroAll(:,idx,neuron_ID),2));
                    
                            
                
                  DFF_idx = DFF_idx + 1 ; 
              else
                  continue
              end
             
              
      end
     
      DFF_mu(jj,:,:) = DFFtemp;
        
     
 end

[df_by_level, df_by_level_timing ] = max(DFF_mu(:,:,:),[],3);

% this function takes a sorted FRA and yields a BF map 