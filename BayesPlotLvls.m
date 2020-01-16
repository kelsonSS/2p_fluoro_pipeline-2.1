function BayesPlotLvls(BayesModelbyLvl)
% this function creates the standard accuracy plot by level for the results 
% of a Bayes Model  - Kelson Shilling-Scrivo 2019 

%input checking 
Test= BayesModelbyLvl ;

if ndims(Test)<4
    error('expected 4d variable: ensure var is in X by Class by Lvl by Rep')
end 
 



% init 


nreps = size(Test,4);
% option 1- saturation of class colors 
Colors(:,:,1) = ...
         [0   0         0   
         0    0.122     0.196;...
         0    0.294     0.482;...
     0.404    0.718     0.922;]
Colors(:,:,2) = ...
   [ 0        0         0;...
     0.412    0.086     0;...
     0.651    0.137     0;...
     1.000    0.545     0.422;]
 Colors(:,:,3) = ...
   [ 0         0         0;...
     0.482    0.337     0;...
     0.710    0.506     0.035;...
     1.000    0.800     0.333;]

% option 2 - parula
% Colors = colormap('parula');
% Colors = Colors([24,34,44,54],:);
% 
% 
% Colors = repmat(Colors,1,1,3);
%      
% % option 3- primary colors  
% Colors(:,:,1) = [  0 0 0;...
%                    1 0 0;...
%                    0 1 0;...
%                    0 0 1];
%
%  Colors = repmat(Colors,1,1,3);                  


%% main 

 Classnames = {'Noise-On','Tone-On','Tone-OFF'};
figure



for class = 1:3 % classes 
    subplot(1,3,class) 
    hold on 
     for lvl = 1:4  % level
        shadedErrorBar([], mean( Test(:,lvl,class,:) , 4),...
             std( Test(:,lvl,class,:) , [] , 4 ) / sqrt(nreps) * 1.96,...
             {'Color', Colors(lvl,:,class)} )    
     end
     
     title(Classnames{class})
     
     if class == 3
     h = gca; 
      legend(h.Children(4:4:16),{'Inf','+20 dB','+10 dB','  0 dB'},'Location','southeast')
     end 
end 
