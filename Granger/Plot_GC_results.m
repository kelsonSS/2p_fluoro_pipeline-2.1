function Plot_GC_results(GC)

% this function takes the GC struct from Collect_GC_results and makes plots
% from it 

% Plot GC fields by level
if isfield(GC,'GC')
    GC = GC.GC;
end 

fn = fieldnames(GC);





    for field = 1:length(fn)
        if ~ strcmp(fn{field}, {'mode','Exptlist','InsufficientCells',...
                                  'GCnumbers_mu','GCnumbers_std', 'GCangles'} )  
    
         data = GC.(fn{field});
         if iscell(data)
             if ~all(cellfun(@kstest,data))
                 warning('GC: %s is not normally distributed!',fn{field}) 
             end 
             
             data = cellfun(@mean,data);
         end 
         
            GCFig(data)
            title([fn{field}]) 
        
             
        
        end 
    end 
    % plot polar plots for GCangles
    data = GC.GCangles;
    lvls = size(data,2);
    rc = numSubplots(lvls);
    figure
    for lvl = 1:lvls
       subplot(rc(1),rc(2), lvl)
        polarhistogram( cat(1,data{:,lvl}), 12)
        title( ['Level ' num2str(lvl)] )
        
    end
        
    
    
    
    
end 

function GCFig(data)
        
        figure
%         bar(mean(data),'k')
%         hold on 
%         errorbar(mean(data), std(data) / sqrt(size(data,1)), 'k.' );
%         
       boxplot(data)

end   
        
        
        
        
    
    
    