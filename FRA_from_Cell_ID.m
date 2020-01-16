function FRAS = FRA_from_Cell_ID(dataDir,IDs,savepath,class)
% takes Cell ids from a directory and plots their FRAs 
%dataDir - list of folders containing the expts
%IDs - Cell of IDs of selected Neurons 
%optional
% savepath - what is the main path to the save location
% class - What subfolder should the files be stored in

assert(length(dataDir) == length(IDs))
if ~exist('class','var')
    class = [] ;
end 
 

for expt = 1:length(dataDir)
    % select traces
    id = IDs{expt};
    F_table = Fluoro_to_Table(dataDir{expt});
    F_table.DFF =  F_table.DFF(:,:,id);
    
    % convert to FRA
    [~,df_by_level] = getTracesAndFRA(F_table);
    
    FRAS{expt} = df_by_level;
    
    freqs = unique(F_table.FreqLevelOrder{:,1});
    levels =  sort(unique(F_table.FreqLevelOrder{:,2}),1,'descend');
    % plotting
    % assuming length(ids) = 20
    
    
    % Average Response Plot
  h1 = figure;
    for nn = 1:length(id)
        try
            subplot(5,4,nn)
            ToneInNoise_MeanTrace(F_table.DFF(:,:,nn),id(nn));
            
        catch
        end
    end
    set(h1, 'Units', 'Inches', 'Position', [0, 0, 11, 8.5],...
             'PaperUnits', 'Inches', 'PaperSize', [11, 8.5])
    
    
    % DFF_Figure 
    h2 = figure;
    for nn = 1:length(id)
        try
            subplot(5,4,nn)
            myFRA(freqs,levels,df_by_level(nn,:), sprintf(' %d ',id(nn) ))
            set(gca,'Ydir','reverse')
        catch
        end
    end
    set(h2, 'Units', 'Inches', 'Position', [0, 0, 11, 8.5],...
             'PaperUnits', 'Inches', 'PaperSize', [11, 8.5])
    
    % saving
    if exist('savepath','var')
        mkdir(fullfile(savepath,'GC_FRA',class));
        mkdir(fullfile(savepath,'Mean_trace',class));
        
        % naming 
        outpath2 = fullfile(savepath,'GC_FRA',class,...
                           sprintf('Experiment_%d.pdf',expt));
                       
         outpath = fullfile(savepath,'Mean_trace',class,...
                           sprintf('Experiment_%d.pdf',expt));
                                     
        
        % saving 
        saveas(h1,outpath)
        
        
        saveas(h2,outpath2)
        close(h1)
        close(h2)
        
        
    end 
    
end
    
    

    