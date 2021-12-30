function resultsAll = getFluoroInfo(TNBehavior,type)
% this function extracts basic info such as the number of active neurons
% and the CellIDS

activeNeurons_name = sprintf('ActiveNeurons_%s',type);
ResponsiveNeurons_name = sprintf('ResponsiveNeurons_%s',type)
for expt_idx = 1:length(TNBehavior)  
    % get psignal info
        TN = TNBehavior{expt_idx}
    
        % munge into correct form for concatenation
        results.n_Neurons = size(TN.DFF,3);
        results.(activeNeurons_name) = sum(TN.active{:,2} >= 1);    
        n_active_levels =  squeeze( sum(...
                                    sum(...
                                        TN.df_by_level_sig{1},...
                                        1),...
                                        2) ) ;
        results.(ResponsiveNeurons_name) = sum( n_active_levels >= 1) ;
        % Package into Output Structure
        if expt_idx == 1
            resultsAll = results;
        else
            resultsAll = ConcatenateStructs(resultsAll,results,1);
        end
        
end
end

