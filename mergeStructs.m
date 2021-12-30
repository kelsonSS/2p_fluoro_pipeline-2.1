function out = mergeStructs(A,B)

fn = fieldnames(B)
for field = 1:length(fn)
    
   A.(fn{field}) = B.(fn{field});
end
    
    
    
    
 out = A;
