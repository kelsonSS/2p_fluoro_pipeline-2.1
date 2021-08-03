function PlotScatterLevels(GroupNames,Levelnames,Title,Style,varargin)
% PlotScatterLevels(GroupNames,Levelnames,Title,varargin)
% creates scatterplots of N datasets at the various levels 
% 
% 
% LevelNames - the names of the levels used 
% GroupNames - the name of the groups used
% style - bar or scatter
% varargin - the data to be plotted, with each group being a seperate entry


figure 
hold on 

% where plotting starts and how much we'll shift each plot by 
increment = .2 ;
start_level = 1;
all_start_levels = []
% to save fig handles 
h = [];

for ii = 1:length(varargin)

    % grab data
    data = varargin{ii};
    n_animals = size(data,2);
    n_levels = size(data,1);
    
    % find where to start plotting levels
    start_level = start_level+increment;
    all_start_levels(ii) = start_level
   levels =  repmat([start_level: n_levels+1]' ,1, n_animals);
   
   % plot scatterplot
   switch Style
       case 'bar'
         h(ii) =bar(levels(:,1),nanmean(data,2)); 
       case 'scatter'
          h(ii) =scatter(levels(:),data(:),'filled');
   end 
   % plot errorbar with 95% CI
   errorbar(levels(:,1), nanmean(data,2), nanstd(data,[],2) / sqrt(n_animals) * 1.96 ,'.');  
   (levels(:,1),nanmean(data,2); 
end

legend(h,GroupNames)
xticks( [ median(all_start_levels) : n_levels+1] )
xticklabels(Levelnames)
xlabel('Levels')

if length(Levelnames) == 1 || ischar(Levelnames)
    xlim( [.5 1.5]+ increment  )
end 

if Title
    title(Title,'interpreter','none')
    saveas(gcf,sprintf('%s.pdf', Title))
end 

end 



