function Stats =  AnalyzeCorrelationsActivePassive(Young_Corr,Old_Corr,SaveName)
    
    if ~exist('SaveName','var')
        SaveName = ''
    end 
    columns = fieldnames(Young_Corr);
    Old_means =[];
    Young_means = [];
    Old_CI = [];
    Young_CI = [];
    Stats = struct();
    
    
    figure 
    hold on 
    suptitle(SaveName)
    rc = numSubplots(length(columns));
    young_above = {};
    young_below = {};
   
    old_above = {};
    old_below = {};
   
    above_fields = {};
    below_fields = {};
    
    above_idx = 1;
    below_idx = 1;
    
    for field = 1:length(columns)
       % extract field
        column_name = columns{field};
        young_field = Young_Corr.(column_name);
        old_field = Old_Corr.(column_name);
        
        % store field
        if contains(column_name,'Above')
            young_above{above_idx} = young_field.gain';
            old_above{above_idx} = old_field.gain';
            above_fields{above_idx} = column_name;
            above_idx = above_idx + 1;
        else 
            young_below{below_idx} = young_field.gain';
            old_below{below_idx} = old_field.gain';
            below_fields{below_idx} = column_name;
            below_idx = below_idx + 1;
        end 
        % plot field 
        subplot(rc(1),rc(2),field)
        plotField(young_field,old_field,column_name)
        
        % perform t-test on field
        Stats.(column_name) = tTestField(young_field,...
                                           old_field);
        % add to arrays for plotting
       Young_means(field) = young_field.mean;
       Young_CI(field) = young_field.CI;
       Old_means(field) = old_field.mean;
       Old_CI(field) = old_field.CI;
                                       
        
    end 
    % save plot 
    plot_name = sprintf('CorrelationsActivePassiveLevels_%s.pdf',SaveName)
    saveas(gcf,plot_name)
    
    
    Stats.r_minus_BF = Compare2AnovaBehavior(young_above,old_above,...
                         {'Young','Old'},...
                         above_fields,...
                         ['CorrelationsActivePassiveLevels_' SaveName '-r_minus_BF']);
                     
    Stats.r_plus_BF =  Compare2AnovaBehavior(young_below,old_below,...
                         {'Young','Old'},...
                         below_fields,...
                         ['CorrelationsActivePassiveLevels_' SaveName '-r_plus_BF']);
    % grouped Analysis
    grouped_means = cat(2,Young_means',Old_means');
    grouped_CI = cat(2,Young_CI',Old_CI');
    
    
    
   %%plotting  
    % plot above and below on seperate plots 
       above_flag = contains(columns,'Above');
       group_names = {'Far-20dB','Far-0dB','Close-20dB','Close-0dB'};
    
  figure
       
    subplot(2,1,1)
   
    PlotGroupedErrorBars(grouped_means(~above_flag,:),grouped_CI(~above_flag,:));
    xticklabels(group_names)
    title([SaveName ' \Deltar+'])
    
       
  
    subplot(2,1,2)
    PlotGroupedErrorBars(grouped_means(above_flag,:),grouped_CI(above_flag,:));
    xticklabels(group_names)
    title([SaveName ' \Deltar-'])
    

    
     plot_name = sprintf('CorrelationsDeltaActivePassiveLevelsBar_%s.pdf',SaveName)
    saveas(gcf,plot_name)
       
     
 
    
    function plotField(young,old,fieldname)
        hold on 
        
       mu = [];
       err = [];
       
       if contains(fieldname,'20')
           line_style = '-'; 
       else 
           line_style = '--';
       end
      %plot young 
       [mu(1),err(1)] = getMeansAndCI(young.passive);
       [mu(2),err(2)] = getMeansAndCI(young.active);      
       
       errorbar(mu,err,['b' line_style])
       % plot old
        mu = [];
       err = [];
       [mu(1),err(1)] = getMeansAndCI(old.passive);
       [mu(2),err(2)] = getMeansAndCI(old.active); 
        
        errorbar(mu,err,['y' line_style])
        
        plot([.9 2.1], [0 0],'k--' )
        xlim([.9 2.1])
        
        xticks([1 2])
        ylim([-.15 .4])
        xticklabels({'Passive','Active'})
        title(fieldname, 'interpreter' ,'None')
    
function [mu CI] = getMeansAndCI(data)
    mu = nanmean(data);
    CI = nanstd(data) / sqrt(length(data)) * 1.96;
            
            
    
    
    function out = tTestField(f1,f2)
       
        % performs a two sample t-test on two distributions following
        % formula on 
        % mathworks.com/help/stats/ttest2.html
        
      n_pooled = min(f1.n,f2.n);
      
      % looking at right sided t-test so want positive difference
        
       mu = abs(f1.mean - f2.mean);
      % get var and std
      
       pooled_var = (f1.var * (f1.n -1) + f2.var * (f2.n-1)) / (f1.n + f2.n - 2);
       
       std = sqrt( pooled_var);
      
       % t statstic and t-test
       t = mu / sqrt(pooled_var/ n_pooled) ;
 
     
       p = tcdf(t, n_pooled,'upper');    
       
       % packaging
       out.t = t;
       out.n = n_pooled;
       out.p = p;
       out.mu = mu;
       out.std = std;
       
       
       
       
        
        
        
        
        
        
    