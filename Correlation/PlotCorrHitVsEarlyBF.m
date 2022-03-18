function out = PlotCorrHitVsEarlyBF(data,SaveName)
 
 
[Hit,Early,field_lists] = extractMungeData(data);

out.Above = Compare2AnovaBehavior(Hit.Above,Early.Above,...
                                 {'Hit','Early'},field_lists.Above,...
                                  [SaveName '-CorrHitVsEarly_dr_minus']);
out.Below = Compare2AnovaBehavior(Hit.Below,Early.Below,...
                                 {'Hit','Early'},field_lists.Below,...
                                  [SaveName '-CorrHitVsEarly_dr_plus']);                              
                              
                              

 
function [Hit_data,Early_data,out_fields] = extractMungeData(data)
   
        
      Hit = data.Hits.Stats;
      Early = data.Early.Stats; 
       
       fields = fieldnames(Hit);
       
       
       above_idx = 1;
       below_idx = 1;
       for field_idx = 1:length(fields)
           
            fn = fields{field_idx};
           
           if contains(fn,'Above')
               Hit_data.Above{above_idx} = Hit.(fn).gain';
               Early_data.Above{above_idx} = Early.(fn).gain';
               out_fields.Above{above_idx} = fn;
               above_idx = above_idx+1;
                      
           elseif contains(fn,'Below')
               Hit_data.Below{below_idx} = Hit.(fn).gain';
               Early_data.Below{below_idx} = Early.(fn).gain';
               out_fields.Below{below_idx} = fn;
               below_idx = below_idx+1;
               
           else
               warning('Encountered field without Above or Below')
           end
               
           
    
       end  
         
        
        
        
        
 
 
 

 
 








