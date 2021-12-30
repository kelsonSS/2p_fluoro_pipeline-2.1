function Out = CorrelationsActive(TNActive)

Out = struct();

varnames = {'Corr','NCorrTotal'}
for expt = 1:length(TNActive)
    expt_corr = Correlations(TNActive{expt});

    for var_idx =1:length(varnames)
        curr_var = varnames{var_idx};
        if expt == 1
        Out.(curr_var) = expt_corr.(curr_var)
        else 
        Out.(curr_var) = cat(1,Out.(curr_var),expt_corr.(curr_var))
        end 
    end 
    


    

end 

