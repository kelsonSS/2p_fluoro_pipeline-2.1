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


if ~exist('activity','var')
    activity = 1; %.05 signifiance
end



try
    handles = In.handles;
catch
    handles = In.Handles;
end
if length(handles)> 1 
    handles = handles(1);
end 

df_by_level = DFFtoFRA(In);

% 
% DFF = In.DFF;
try
 active_list = In.Active.Activity;
catch 
   active_list = In.active.Activity; 
end 
 FreqLevels = unique(In.FreqLevelOrder);
 Freqs =In.FreqLevelOrder{:,1};
 Levels= In.FreqLevelOrder{:,2};
 uLevels= unique(Levels);
 uFreqs = unique(Freqs);
 
 soundon = handles.PreStimSilence * handles.pfs;
 soundoff = soundon + handles.PrimaryDuration + handles.pfs; 


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
         n_fig = n_fig+1; 
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
    



    
    