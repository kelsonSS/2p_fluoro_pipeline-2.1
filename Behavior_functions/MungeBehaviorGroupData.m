function Out = MungeBehaviorGroupData(behavior)
%  function to group all animals behavior for an entire group for downstream
% statistical analysis 

% behavior- is a struct with each entry containing named entries as outputs
% from BehavioralAnalysisByAnimal
%  ie Behavior.animal_01.resultsAll 
% animal name is not critical but resultsAll and its relative location is 


SNRs = [20 10 0 ];  
AnimalIDs  =  fieldnames(behavior);



    for expt_idx = 1:length(AnimalIDs)
        Animal_ID = AnimalIDs{expt_idx};
        if strcmp(Animal_ID,'group')
            continue
        end 
% grab expt data for current animal
    expt = behavior.(Animal_ID).resultsAll.Combined;
    
% Subset Data to only those in correct SNRs 
    expt_subset = SubsetExpt(expt,SNRs);

 % Package into Output Structure  	
    if expt_idx == 1
        Out = expt_subset;
    else
        Out = ConcatenateStructs(Out,expt_subset)
    end 
    
    end 
   Out.uSNRs = SNRs';

end 
    

function  Out = SubsetExpt(expt,SNRs)

Out = struct();
% get indicies
for curr_SNR = 1:length(SNRs)
    
  idx = find(expt.SNR' == SNRs(curr_SNR));

 
  Out = fillOutWithValue(Out,expt,idx);

end 
end 

function MainStruct =  ConcatenateStructs(MainStruct,SecondStruct)
% given two structures with the same fields, will attempt to concanenate
% each field columnwise 

fields = fieldnames(MainStruct);

for f_idx = 1:length(fields)
    % determine if fields exist in Output
       c_field = fields{f_idx};
       MainStruct.(c_field) = cat(2,MainStruct.(c_field),...
                                   SecondStruct.(c_field) );

end

end 
    

function Out= fillOutWithValue(Out,expt,idx)


fields = fieldnames(expt);
Out_fields = fieldnames(Out);

for f_idx = 1:length(fields)

    c_field = fields{f_idx}; 
    
    % extract field value at idx with Nan if idx doesn't exist
    if idx
        field_value = expt.(c_field)(idx);
    else 
        field_value = nan;
    end 
    
    % determine if fields exist in Output
    field_in_output_flg = any(strcmp(c_field,Out_fields));
    % add or concatenate values to output structure
    if field_in_output_flg
       Out.(c_field) = cat(1,Out.(c_field),field_value);
    else
       Out.(c_field) = field_value;
    end 
end



end


