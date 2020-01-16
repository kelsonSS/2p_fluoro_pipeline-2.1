function PlotFRAS(DFF,activity,save)

% This function takes matrix in the form of a m x n Neuron X condition matrix 
% with the values indicating the mean fluorescence for that neuron and condition

% this function uses the associated myFRA and numSubplots function to plot multiple FRAs
% onto the same axis 

% outpath1 = 'C:\Users\Kelson\Dropbox\Passive_paper\PV\FRAs\FRA'
% outpath2 = 'C:\Users\Kelson\Dropbox\Passive_paper\PV\FRAs\MeanTrace'

handles = DFF.handles;
active_list = DFF.active.Activity;
Tones = unique(DFF.FreqLevelOrder{:,1});
Levels = unique(DFF.FreqLevelOrder{:,2});
soundon = handles.PreStimSilence * handles.pfs;
soundoff = soundon + handles.PrimaryDuration + handles.pfs; 

    if ~exist('activity','var')
    activity = 1; %.05 signifiance 
    end 

df_by_level = mean(DFF.DFF)

n = 20; % number of plots on same figure
%nTones = 7;
%nLevels =  4;
rc = numSubplots(n);

 h1 = figure
 h2 = figure
nn = 0 
n_fig = 1 
for ii =1:size(df_by_level)
    
    if active_list(ii) >= activity
      nn = nn+1 ;
   figure(h1)
      subplot(rc(2),rc(1),mod(nn-1,n)+1)
    
     myFRA(Tones,Levels,df_by_level(ii,:),ii)
    figure(h2)
     subplot(rc(2),rc(1),mod(nn-1,n)+1)
     hold on 
  
     plot(squeeze(nanmean(DFF.DFF(:,:,ii),2)));
     
     xlm = get(gca,'Xlim');
     plot([xlm(1), xlm(2)],[0, 0])
     hold off
    end 
    
     if mod(nn,n) == 0
      SaveStyleFig(h1,outpath1,n_fig)    
      SaveStyleFig(h1,outpath2,n_fig)   
         h1 = figure;
         h2 = figure;
         nn= 0;
         n_fig = n_fig+1 
     end 
     
  
end 
% for last figure after loop
  SaveStyleFig(h1,outpath1,n_fig,save)    
  SaveStyleFig(h2,outpath1,n_fig,save) 
end 

function SaveStyleFig(h1,outpath,n_fig,save)

set(h1, 'Units', 'Inches', 'Position', [0, 0, 11, 8.5],...
             'PaperUnits', 'Inches', 'PaperSize', [11, 8.5])
if save
    saveas(h1,fullfile(outpath,sprintf('figure%d.pdf',n_fig)))         
end 
end 
    



    
    