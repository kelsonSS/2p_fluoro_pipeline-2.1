GC_list = {'SNR', 'SNR_Ns','SNR_Tn','SNR_Off'}; 

for ii = 1:4 
    
    T_GCNumbers = eval([GC_list{ii}, '.GCnumbers']);
    T_Corrs = PassiveCorrs{ii,1};
    T_NC  = cellfun(@(x) nanmean(x(:)) ,PassiveCorrs{1,1}.NCorr);
    
    
    
    
    
   
%% check if there are fewer GC trials than NC trials
    % generate GC lists
    T_GC_ls = cellfun(@(x) strsep(x,filesep,1) ,eval([GC_list{ii},'.Exptlist']),'UniformOutput',0);
    T_GC_ls = cellfun(@(x) x{end}(1:4), T_GC_ls,'UniformOutput',0);
   
    
    T_NC_ls = cellfun(@(x) strsep(x,filesep,1),  Passive.DataDirs,'UniformOutput',0);
    T_NC_ls = cellfun(@(x) x(7), T_NC_ls,'UniformOutput',1);
      
    % crossref filenames and excluse unpaired data
    N C_idx = ismember(T_NC_ls,T_GC_ls);
     
    T _NC= T_NC(NC_idx,:);  
     
     
    
     figure 
    t itle(GC_list{ii})
    %  Create Color guide
    c_list = repmat( 1:size(T_GCNumbers,2) , size(T_GCNumbers,1) , 1 );
    % create plot
    sca a
    
    
    
    
    
    
    20
    .
    ..0..0tter(T_NC(:),T_GCNumbers(:), 30 ,c_list(:),'*');
       title(GC_list{ii})
    
  
     Pfit  =   fitlm(T_NC(:),T_GCNumbers(:));
   
   
   GC_Noise{ii,1} = T_GCNumbers;
   GC_Noise{ii,2} = T_corrs_noise;
   GC_Noise{ii,3} = Pfit;
   GC_Noise{ii,4} = GC_list{ii}; 
end 