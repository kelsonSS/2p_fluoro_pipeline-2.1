function PlotFRAS(In,activity,Save)

% This function takes matrix in the form of a m x n Neuron X condition matrix 
% with the values indicating the mean fluorescence for that neuron and condition

% this function uses the associated myFRA and numSubplots function to plot multiple FRAs
% onto the same axis 

 outpath1 = '\\Vault3\Data\Kelson\Aging\FRAs\FRA';
 outpath2 = '\\Vault3\Data\Kelson\Aging\FRAs\MeanTrace';
 
% mkdir(outpath1)
% mkdir(outpath2)

 if ~exist('Save','var')
     Save = 0;
 end 

try
handles = In.handles;
catch 
 handles = In.Handles;
end 
if length(handles)> 1 
    handles = handles(1);
end 


DFF = In.DFF;
active_list = In.Active.Activity;
FreqLevels = unique(In.FreqLevelOrder);
Freqs =In.FreqLevelOrder{:,1};
Levels= In.FreqLevelOrder{:,2};
uLevels= unique(Levels);
uFreqs = unique(Freqs);

soundon = handles.PreStimSilence * handles.pfs;
soundoff = soundon + handles.PrimaryDuration + handles.pfs; 

    if ~exist('activity','var')
    activity = 1; %.05 signifiance 
    end 

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
    

n = 20; % number of plots on same figure
%nTones = 7;
%nLevels =  4;
rc = numSubplots(n);

 h1 = figure;
 h2 = figure;
nn = 0 ;
n_fig = 1; 
for ii =1:size(df_by_level)
    
    if active_list(ii) >= activity
      nn = nn+1 ;
   figure(h1)
      subplot(rc(2),rc(1),mod(nn-1,n)+1)
    
     myFRA(uFreqs,uLevels,df_by_level(ii,:),ii)
    figure(h2)
     subplot(rc(2),rc(1),mod(nn-1,n)+1)
     hold on 
  
     plot(squeeze(nanmean(In.DFF(:,:,ii),2)));
     
     xlm = get(gca,'Xlim');
     plot([xlm(1), xlm(2)],[0, 0])
     hold off
    
      if mod(nn,n) == 0 
         if Save
             SaveStyleFig(h1,outpath1,n_fig,Save)
             SaveStyleFig(h1,outpath2,n_fig,Save)
         end
         h1 = figure;
         h2 = figure;
         nn= 0;
         n_fig = n_fig+1 
     end 
    
    
    
    end 

    
     
  
end 
% for last figure after loop
  SaveStyleFig(h1,outpath1,n_fig,Save)    
  SaveStyleFig(h2,outpath1,n_fig,Save) 
end 

function SaveStyleFig(h1,outpath,n_fig,Save)

set(h1, 'Units', 'Inches', 'Position', [0, 0, 11, 8.5],...
             'PaperUnits', 'Inches', 'PaperSize', [11, 8.5])
if Save
    saveas(h1,fullfile(outpath,sprintf('figure%d.pdf',n_fig)))         
end 
end 
    



    
    