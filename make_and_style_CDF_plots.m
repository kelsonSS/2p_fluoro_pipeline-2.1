% make_and_style_CDF_plots
% requires to runn the Corr2D script 

figure
hold on
cdfplot(Corr2D_Flat{1,1})
cdfplot(Corr2D_Flat{2,2})
cdfplot(Corr2D_Flat{3,3})
cdfstyle(gcf,'Noise Correlations',{'Noise',...
                                   'Tone',...
                                   'Offset'} );
    
figure
hold on 
cdfplot(Corr2D_Flat{1,1})
cdfplot(Corr2D_Flat{1,2})
cdfplot(Corr2D_Flat{1,3})
cdfstyle(gcf,'Noise Correlations',{'Noise:Noise',...
                                   'Noise:Tone',...
                                   'Noise:Offset'} );


figure
hold on 
cdfplot(Corr2D_Flat{2,1})
cdfplot(Corr2D_Flat{2,2})
cdfplot(Corr2D_Flat{2,3})
cdfstyle(gcf,'Noise Correlations',{'Tone:Noise',...
                                   'Tone:Tone',...
                                   'Tone:Offset'} );
figure
hold on 
cdfplot(Corr2D_Flat{3,1})
cdfplot(Corr2D_Flat{3,2})
cdfplot(Corr2D_Flat{3,3})
cdfstyle(gcf,'Noise Correlations',{'Offset:Noise',...
                                   'Offset:Tone',...
                                   'Offset:Offset'} );

figure
Class_Corrs = {Corr2D_Flat{1,1},Corr2D_Flat{2,2},Corr2D_Flat{3,3}};
[temp.p,temp.mu,temp.stats] = anovan(Flatten_Corr(Class_Corrs),...
    Flatten_Corr(Corr2D_class_ID));
temp.sig = multcompare(temp.stats);


Class_Corr_Stats = temp;
clear temp 

                               
                               
                               
 function cdfstyle(handle,fig_title,legend_list) 
 
 h = handle; 
 figure(h)
 title(fig_title)
 legend(legend_list)
 axis tight 
 ylabel('Cumulative percentage(%)')
 xlabel('Correlation (A.U)')
 
 
 end 
