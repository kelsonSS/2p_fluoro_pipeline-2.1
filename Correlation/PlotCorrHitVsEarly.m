function out = PlotCorrHitVsEarly(data,SaveName)
 
 
[Hit,Early,field_lists] = extractMungeData(data);

out.dr_minus = Compare2AnovaBehavior(Hit.Above,Early.Above,...
                                 {'Hit','Early'},field_lists,...
                                  [SaveName '-CorrHitVsEarly_dr_minus']);
out.dr_plus = Compare2AnovaBehavior(Hit.Below,Early.Below,...
                                 {'Hit','Early'},field_lists,...
                                  [SaveName '-CorrHitVsEarly_dr_plus']);                              
                              
                              

 
function [Hit_data,Early_data,out_fields] = extractMungeData(data)
   
      out_fields = {'20db','0dB'};
      Hit = data.Hits.Stats;
      Early = data.Early.Stats; 
      Hit_data = struct();
      Early_data = struct();
      
       fields = fieldnames(Hit);
       
       
       
       above_idx = 1;
       below_idx = 1;
       for field_idx = 1:length(fields)
           
            fn = fields{field_idx};
           
           if contains(fn,'Above')
               out_field = 'Above';
           else
               out_field = 'Below';
           end 
            
           if contains(fn,'20')
               cell_idx = 1;
           else
               cell_idx = 2;
           end
           
           
           Hit_data = addData(Hit_data,out_field,cell_idx,Hit.(fn).gain');
               Early_data = addData(Early_data,out_field,cell_idx,Early.(fn).gain');
           
              
          
               
           
    
       end  
         
        
        
        
        
 
 
    function  out = addData(out,out_field,cell_idx,data_to_add)
        
        if ~isfield(out,out_field)
            out.(out_field) = cell(1,2);
            
        end   
        
        old_data = out.(out_field){cell_idx};
        
        
        out.(out_field){cell_idx} = cat(2,old_data,data_to_add);
        
        
            
 

 
 








