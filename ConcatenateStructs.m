
function MainStruct =  ConcatenateStructs(MainStruct,SecondStruct,dim)
% given two structures with the same fields, will attempt to concanenate
% each field columnwise 
if ~exist('dim','var')
    dim = 2;
end

fields = fieldnames(MainStruct);

for f_idx = 1:length(fields)
    % determine if fields exist in Output
       c_field = fields{f_idx};
       MainStruct.(c_field) = cat(dim,MainStruct.(c_field),...
                                   SecondStruct.(c_field) );

end

end 