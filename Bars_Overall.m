function Bars_Overall(mu,CI,name,sig)
       % takes a mean accross levels for both the correlation and
       % confidence interval then plots them

all_mu = mean(mu);
all_CI = mean(CI);
% plots average across levels

figure
errorbar(all_mu,all_CI, '.')
hold on
bar(all_mu,'k')
set(gca,'ylim',[0 .3]) 
if exist('sig','var')
    
    sig_idx = sig(:,end)<=.05 ;
    if any(sig_idx)
        sigstar(sig(sig_idx,1:2),sig(sig_idx,end),[.05,.01,.001]);
    end
end 


hold off
title(sprintf('%s Correlations',name),'Interpreter','none')
set(gca,'XTick',[1])
set(gca,'XTickLabel',{'Average'})
