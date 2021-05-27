 function Bars_by_level(mu,CI,name,sig)
       % takes a mean and 95% confidence interval and plots them 
        % plotting- Ncorr Bars
figure
errorbar(mu,CI, '.')
hold on
bar(mu,'k')
set(gca,'ylim',[0 .5]) 
if exist('sig','var')
    
    sig_idx = sig(:,end)<=.05 ;
    if any(sig_idx)
        sigstar(sig(sig_idx,1:2),sig(sig_idx,end),[.05,.01,.001]);
    end
end 


hold off
title(sprintf('%s Correlations by level',name),'Interpreter','none')
set(gca,'XTick',[1:3])
set(gca,'XTickLabel',{'+20','+10','0'})



