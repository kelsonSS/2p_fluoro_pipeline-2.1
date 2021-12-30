function myFRA(tones,levels,fra,nn,ax)
% this function plots the Frequency response Area of an auditory neuron to
% tone stimuli presented at different levels
% 
if ~exist('ax','var')
    ax = true;
end 
  
 t= length(tones);
 L= length(levels);

 if t ==  1 && L == 1
     t = tones;
     L = levels;
 end 
 
 % normalize  df responses to 1 
fra ./ max( fra(:) )  ;
 
%% reshape into FRA, add zeros to bottom row since pcolor colors by top
 % left position 
 if isrow(fra) || iscolumn(fra)
   fra = reshape(fra,L,t);
 end 
 
 if size(fra,2) == L && size(fra,1) == t
     fra = permute(fra,[2,1]);
 end 
 
  fra = [[fra ; nan(1,t)], nan(L+1,1)];
  
  
  
  
  
  %% create figure 
   h = pcolor(fra) 
   set(h, 'EdgeColor','none')
   set(gca,'Ydir','Reverse')    

   
 %% Format Figure 
 
 % axis ij; %levels = levels(end:-1:1);
 if ax
     set(gca,'XTick',[1:t])
     set(gca,'YTick',[1:L])
     set(gca,'XTickLabel',round(tones/1000,0))
     set(gca,'YTickLabel',levels)
     title(['Neuron ' num2str(nn)])
     yticks(yticks + .5)
     xticks(xticks +.5)
     xlabel('Frequency')
     ylabel('Level')
     colorbar
 else
     axis off
 end
     
             

        
        
        

