function freqs_prop = PlotBFDistribution(DF,lvl)

    if ~exist('lvl','var')
        lvl = 0
    end 
  df_by_level = [];
  df_by_level_sig = [];

            for ii = 1:length(DF.df_by_level)
                try
                df_by_level = cat(3,df_by_level,DF.df_by_level{ii});
                df_by_level_sig = cat(3,df_by_level_sig,DF.df_by_level_sig{ii});
                catch
                    continue
                end 
            end
            
   %df_by_level = df_by_level .* df_by_level_sig;
   
   if size(df_by_level,1) > 1 % if there are multiple levels 
   
       if lvl
           df_by_level = squeeze(df_by_level(lvl,:,:));
       else
           
           df_by_level = squeeze(max(df_by_level));
       end
       
   else
       df_by_level  = squeeze(df_by_level);
   end
   
   n_freqs = size(df_by_level,1);
   n_neurons = size(df_by_level,2);
   
   [~,freqs] = max(df_by_level);
   freqs_count = histcounts(freqs,n_freqs);
   freqs_prop = freqs_count ./ n_neurons;
   
   
   figure;bar(freqs_prop)
   
   