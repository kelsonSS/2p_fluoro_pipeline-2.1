function FullAnimalInfo = CreateFullAnimalExperimentSheet(TNBehavior)

% first we create a individual table of the base animal info
AnimalInfo = getAnimalInfo(TNBehavior(:,4));
% grab behavioral info
BehavioralInfo = struct2table(BehavioralAnalysisByExperiment(TNBehavior(:,4)));
BehavioralLevelInfo = getLevelData(BehavioralInfo, [20,10,0]);

% grab imaging info
PassiveInfo = struct2table(getFluoroInfo(TNBehavior(:,1),'Passive'));
ActiveInfo  = struct2table(getFluoroInfo(TNBehavior(:,2),'Active'));



DetectionInfo = DetectionAnalysis(TNBehavior(:,2));
DetectionInfo = structfun(@(x) x', DetectionInfo,'UniformOutput',0);
DetectionInfo = struct2table(DetectionInfo);
DetectionInfo = DetectionInfo(:,1:3);
LickResponseTimes = LoadLickingBehavior(TNBehavior(:,4));
LickResponseTimes = table(LickResponseTimes);

BF_Field = table(BestFrequencyAnalysisAll(TNBehavior(:,1))',...
                 'VariableNames',{'BF_Field'});
% active and passive have the same number of neurons so we need to drop
% duplicated row 
ActiveInfo.n_Neurons = [];

% create the full table of all subtables  
FullAnimalInfo = [ AnimalInfo , BehavioralInfo,BehavioralLevelInfo,...
                PassiveInfo,BF_Field , ActiveInfo,DetectionInfo,LickResponseTimes];

            
            
    function out= getLevelData(data, SNRs_to_use)
      
        out = {};
    vars_to_use =   find(contains(data.Properties.VariableNames,'Level'));     
    SNRs = data.SNRs ;
    
    if isempty(vars_to_use)
        return
    end 
    for var_idx = 1:length(vars_to_use)
      
            curr_var = vars_to_use(var_idx);
          temp_table = createTableWithLevels(data(:,curr_var), SNRs,SNRs_to_use);
          out = [out, temp_table];

    end 
    
   
        function temp_table = createTableWithLevels(curr_table,SNRs,SNRs_to_use)
            
           baseName = curr_table.Properties.VariableNames;
           VarNamesLevels =  arrayfun(@(x) [baseName{1}, '_' , num2str(x)],...
                             SNRs_to_use,'UniformOutput',0);
           
          n_expts = height(curr_table)
          temp_table = table()
          
          for lvl_idx = 1:length(SNRs_to_use)
              
              current_SNR = SNRs_to_use(lvl_idx);
              % preallocate
              temp_table.(VarNamesLevels{lvl_idx}) = nan(n_expts,1);
              
              % fill
              for expt_idx = 1:n_expts
                  
                expt_SNR_idx =   find( SNRs{expt_idx} ==current_SNR);
                
                if isempty(expt_SNR_idx)
                    continue
                else 
                    val_expt = curr_table{expt_idx,baseName}{1}(expt_SNR_idx);
                    temp_table{expt_idx,lvl_idx} = val_expt;
                end 
              
              
              end 
          end
          
          
           
           
           
           
           
            
    
        
        
        
    
  
     