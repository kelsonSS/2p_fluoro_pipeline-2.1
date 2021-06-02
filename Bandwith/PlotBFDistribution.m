function freqs_prop = PlotBFDistribution(DF,lvl,by_expt)

    if ~exist('lvl','var')
        lvl = 0;
    end 
    
    if ~exist('by_expt','var')
        by_expt = 0 ;
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

%   only look at cells with significant responses 
    
    %df_by_level = df_by_level .* df_by_level_sig
     df_by_level = df_by_level(:,:,DF.active{:,2}>0) ;
     expt_list = DF.experiment_list(DF.active{:,2}>0);
     
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
  
   if by_expt
       expt_id =1;
       for expt = 1:max(DF.experiment_list)
          
           
           df_expt = df_by_level(:,expt_list == expt);
           n_neurons = size(df_expt,2);
           if isempty(df_expt)
               continue
           end
           [ ~,freqs] = max(df_expt);
           freqs_count = histcounts(freqs,n_freqs);
           freqs_prop_t = freqs_count ./ n_neurons;
           
           freqs_prop(:,expt_id) = freqs_prop_t;
           expt_id = expt_id +1
           clear freqs_prop_t
       end 
           
           
           figure;errorbar([],mean(freqs_prop'),std(freqs_prop')/sqrt(size(freqs_prop,2)) * 1.93,'.')
           hold on
           bar(mean(freqs_prop') , 'k' )
           
           figure;bar(freqs_prop)
           
       
   else 
    
       
   
   [~,freqs] = max(df_by_level);
   freqs_count = histcounts(freqs,n_freqs);
   freqs_prop = freqs_count ./ n_neurons;
   
   
   figure;bar(freqs_prop)
   
        end 