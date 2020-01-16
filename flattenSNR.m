function SNR_flat = flattenSNR(SNR)

flds = fieldnames(SNR);
for name_idx = 1:length(flds)
    name =  flds( name_idx);
    name = name{1};
    [r,c] = size( SNR.( name ));
    
    if r >1 && c > 1 && iscell( SNR. (name) )
        temp = SNR.( name );
       for cc = 1:c 
           temp_col = [];
            for rr = 1:r
                temp_col = cat(1,temp_col,temp{rr,cc} ) ;
            end
               out{cc} = temp_col;
       end 
        
        
        SNR_flat.(name) = out;
    else 
        SNR_flat.( name ) = SNR.( name );
    end 
end 