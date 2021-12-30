
function  Out = SubsetExpt(expt,SNRs)

Out = struct();
% get indicies
for curr_SNR = 1:length(SNRs)
    
  idx = find(expt.SNR' == SNRs(curr_SNR));

 
  Out = fillOutWithValue(Out,expt,idx);

end 
end 



function Out= fillOutWithValue(Out,expt,idx)


fields = fieldnames(expt);
Out_fields = fieldnames(Out);

for f_idx = 1:length(fields)

    c_field = fields{f_idx}; 
    
    % extract field value at idx with Nan if idx doesn't exist
    if  ~idx
        field_value = nan;
    % if there is only 1 element, do not try to subset
    elseif numel(expt.(c_field)) == 1 
          field_value = expt.(c_field);
    else 
     % subset to only desired levels 
        field_value = expt.(c_field)(idx);
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
