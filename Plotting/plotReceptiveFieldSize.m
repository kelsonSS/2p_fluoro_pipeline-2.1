function rf_sum = plotReceptiveFieldSize(Passive)

   % Creates receptive field area maps of responsive neurons across levels
  
   % if there are no classes make a class of all active neurons 
   if ~ isfield(Passive,'Classes')
       Passive.Classes =  Passive.active{:,2}> 0 ;
   end
   
   % intialize info
   classes = unique(Passive.Classes);
   L = length(unique(Passive.FreqLevelOrder.Levels));
   F = length(unique(Passive.FreqLevelOrder.Freqs));
  
   % FRAs are in Cell x condition order  
   FRAs = DFFtoFRA(Passive);
   % normalize each FRA such that the maximum response to its preferred
   % stimulus is 1 
   FRAs = FRAs./max(FRAs,[],2);
   
   
   figure;
   hold on 
  for class = :max(classes) % don't plot nonactive neurons which are class 0  
      
      % subset FRAS by class
      FRA = FRAs(Passive.Classes == class,:);
      
      % reshape to Level by Frequency matrix      
      FRA = reshape(FRA,size(FRA,1),L,F);
      
      RF = squeeze(sum(FRA,3));
      
      
      
      errorbar(mean(RF),std(RF))
  end    
      % Calculate Receptive field size 
      
      
      
      
      
      
  end 