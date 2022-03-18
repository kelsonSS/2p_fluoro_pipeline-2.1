function  out = FullCompareCorrs(passive,active,close_idx,PlotColor,SaveName,Template)

if ~exist('Template','var')
    Template = []
end

[~,passive_all_mu] = MungeCorr2Mat(passive.NCorrAll);
[~, active_all_mu] = MungeCorr2Mat(active.NCorrAll);

if ~isempty(Template)
    out.All = NikCompareCorrs(passive_all_mu,active_all_mu,close_idx,PlotColor,Template.All);
else
    out.All = NikCompareCorrs(passive_all_mu,active_all_mu,close_idx,PlotColor,[]);
end 

% figure;
% hold on
% cdfplot(passive_all2)
% cdfplot(active_all2)
% legend({'Passive','Active'})
% title('Noise Correlations')
if SaveName
saveas(gcf, sprintf('%s_NoiseCorrelations_Nik.pdf',SaveName))
end 

%%by Level 
LevelNames = {'SNR20','SNR10','SNR0'};

for lvl =1:length(LevelNames)
 [~,active_lvl_mu] = MungeCorr2Mat(active.NCorr(:,lvl));
 
 LevelName = LevelNames{lvl};  
 if ~isempty(Template) 
 out.(LevelName) =  NikCompareCorrs(passive_all_mu,active_lvl_mu,close_idx,PlotColor,Template.(LevelName));
 
 else
 out.(LevelName) =  NikCompareCorrs(passive_all_mu,active_lvl_mu,close_idx,PlotColor,[]);
 Template.SNR10 = out.(LevelName);
 Template.SNR0 = out.(LevelName);
 end
 if SaveName
saveas(gcf, sprintf('%s_NoiseCorrelations_Nik_%s.pdf',SaveName,LevelName))
 end 
end 

plotLevels(out,LevelNames)
 if SaveName
saveas(gcf, sprintf('CorrelationActivePassive_%s_Levels.pdf',SaveName))
 end 


[Close_20_Above] =getSubtractedStats(...
                          out.SNR20.BF.Close.Passive.Above,...
                          out.SNR20.BF.Close.Active.Above);
[Close_20_Below] =getSubtractedStats(...
                          out.SNR20.BF.Close.Passive.Below,...
                          out.SNR20.BF.Close.Active.Below);

 [Close_0_Above] =getSubtractedStats(...
                          out.SNR0.BF.Close.Passive.Above,...
                          out.SNR0.BF.Close.Active.Above)  ;
                      
  [Close_0_Below] =getSubtractedStats(...
                          out.SNR0.BF.Close.Passive.Below,...
                          out.SNR0.BF.Close.Active.Below)   ;                   

 Close_means = [Close_20_Above.mean,Close_0_Above.mean,...
                Close_20_Below.mean,Close_0_Below.mean]';
 
 Close_CIs = [Close_20_Above.CI,Close_0_Above.CI,...
                Close_20_Below.CI,Close_0_Below.CI]';
            
  
% repeat for far 
[Far_20_Above] =getSubtractedStats(...
                          out.SNR20.BF.Far.Passive.Above,...
                          out.SNR20.BF.Far.Active.Above);
[Far_20_Below] =getSubtractedStats(...
                          out.SNR20.BF.Far.Passive.Below,...
                          out.SNR20.BF.Far.Active.Below);

 [Far_0_Above] =getSubtractedStats(...
                          out.SNR0.BF.Far.Passive.Above,...
                          out.SNR0.BF.Far.Active.Above)  ;
                      
  [Far_0_Below] =getSubtractedStats(...
                          out.SNR0.BF.Far.Passive.Below,...
                          out.SNR0.BF.Far.Active.Below)   ;                   

 Far_means = [Far_20_Above.mean,Far_0_Above.mean,...
                Far_20_Below.mean,Far_0_Below.mean]';
 
 Far_CIs = [Far_20_Above.CI,Far_0_Above.CI,...
                Far_20_Below.CI,Far_0_Below.CI]';

            
Means = cat(2,Close_means,Far_means);
CIs = cat(2,Close_CIs,Far_CIs);



figure
PlotGroupedErrorBars(Means,CIs)
xticklabels({' \Deltar- 20dB',' \Deltar- 0dB', ' \Deltar+ 20dB',' \Deltar+ 0B'})
legend({'BF <.5', 'BF >.5'}) 
ylim([-.3 .2])
if SaveName
    saveas(gcf, sprintf('%s_NoiseCorrelations_Nik_Bar_%s.pdf',SaveName,'All'))
end 
      
Stats.Far_20_Above = Far_20_Above;
Stats.Far_0_Above = Far_0_Above;
Stats.Far_20_Below = Far_20_Below;
Stats.Far_0_Below = Far_0_Below;

Stats.Close_20_Above = Close_20_Above;
Stats.Close_0_Above = Close_0_Above;
Stats.Close_20_Below = Close_20_Below;
Stats.Close_0_Below = Close_0_Below;

out.Stats = Stats;    
    
    
    
    
      

function plotLevels(D,LevelNames)
  figure
  n_subplots = length(LevelNames);
  rc = numSubplots(n_subplots);

  for lvl_idx = 1:n_subplots
      
      lvl_title = LevelNames{lvl_idx};
      
     subplot(rc(1),rc(2),lvl_idx)
     
      plotLevel(D.( lvl_title))
      hold on 
      title(lvl_title)
            
            
end 
       
    function plotLevel(d)
        
        
        passive = d.Passive;
        active = d.Active;
        r_plus_idx = d.r_plus_cells;
        r_minus_idx = d.r_minus_cells;
        
        passive_plus = passive(r_plus_idx);
        passive_minus = passive(r_minus_idx);
        
        active_plus = active(r_plus_idx);
        active_minus = active(r_minus_idx);
        
       
       data= {passive_plus,passive_minus;active_plus,active_minus};
       
       means = cellfun(@nanmean, data);
       CIs = cellfun(@get95CI,data);
      
       errorbar([1.3,1.3;2.3,2.3],means,CIs,'LineWidth',2)
       hold on 
       
       plot([1.3,2.3],[.1,.1], 'k--')
       
       ylim([-.15 .3])
       xticks([1.3,2.3])
       xticklabels({'Passive','Active'})
       legend({'r_minus','r_plus'})
       
       
        
        

        function CI = get95CI(data)
            
            CI = nanstd(data) / sqrt(length(data)) * 1.96 ;
            

function out = subsetActiveCorrs(Reference,Corrs,SNRs,uLevels_expt)

    len_SNRs=  length(SNRs);
    len_corrs = length(Reference{1}(:));
    out = nan(len_corrs,len_SNRs);
    
    
for ii = 1:length(SNRs)
    
    idx = find( uLevels_expt == SNRs(ii));
    
    if (~isempty(idx)) && length(Corrs)>=idx && ( length(Corrs{idx}(:) ) == len_corrs )
        out(:,ii) = Corrs{idx}(:);
    end 
end 
        
    

function [m,m_mu] = MungeCorr2Mat(corrs)
 m = cell2mat(cellfun(@(x) x(:),corrs,'UniformOutput',0));
 m_mu = cell2mat(cellfun(@(x) nanmean(x)', corrs,'UniformOutput',0));



    function  NoiseCorrResults = NikCompareCorrs(rPDFF, rHDFF,close_idx,PlotColor,All)

 if ~exist('All','var')
     All = [];
 end 
NoiseCorrResults = [];
NoiseCorrResults.Passive = rPDFF;
NoiseCorrResults.Active  = rHDFF;

Fig = figure

subplot(2,2,1)

n = min(length(rHDFF),length(rPDFF));

rPDFF = rPDFF(1:n);
rHDFF = rHDFF(1:n);

bins = [-.6:.05:.6];

%bins = linspace(min(rPDFF),max(rPDFF),20);

%Passive
CorrHist = hist(rPDFF,bins);
CorrHist = CorrHist./sum(CorrHist);
plot(bins,cumsum(CorrHist),'color',[.33 .33 .33],'linewidth',2)
hold on
%Hit
%bins = linspace(min(rHDFF),max(rHDFF),20);

CorrHist = hist(rHDFF,bins);
CorrHist = CorrHist./sum(CorrHist);
plot(bins,cumsum(CorrHist),'b','linewidth',2)
xlim([-.6 .6])
%plot([0 0],[0 1],'k')
xlabel('Noise Corr.')
ylabel('Probability')
set(gca,'fontsize',8)
title('Noise Correlation (r) CDF')

% %%%%%%%% Statisitical tests for hit vs pass %%%%%%%% 
% %hit
% A = rHDFF;
% [h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
% NoiseCorrResults.Hit.mu=nanmean(A(~isnan(A)));
% NoiseCorrResults.Hit.std=2*std(A(~isnan(A)))./sqrt(length((~isnan(A))));
% NoiseCorrResults.Hit.bootsig001=[h p];
% %passive
% A = rPDFF;
% [h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
% NoiseCorrResults.Pass.mu=nanmean(A(~isnan(A)));
% NoiseCorrResults.Pass.std=2*std(A(~isnan(A)))./sqrt(length((~isnan(A))));
% NoiseCorrResults.Pass.bootsig001=[h p];
% %Hit vs passive
% A = rHDFF;
% B = rPDFF;
% [h p]=bootsig(A(~isnan(A)), B(~isnan(B)), 10000, 0.001);
% NoiseCorrResults.HitPass.mu=nanmean(A(~isnan(A)))-nanmean(B(~isnan(B)));
% NoiseCorrResults.HitPass.std=2*nanstd(A-B)./sqrt(length(A-B));
% NoiseCorrResults.HitPass.bootsig=[h p];

%%%%%%%% Find Hit-Passive CDF diff inflection point
CorrHist = hist(rPDFF,bins);
PCorrHist = cumsum(CorrHist./sum(CorrHist));
CorrHist = hist(rHDFF,bins);
HCorrHist = cumsum(CorrHist./sum(CorrHist));
diffCorrHist = HCorrHist-PCorrHist;
[junk m]=max(diff(diffCorrHist));
InflP = bins(m);
NoiseCorrResults.InflP = InflP;

if ~isempty(All)
    KeptIdxAbove = All.r_minus_cells;
    KeptIdxBelow = All.r_plus_cells;
else
    KeptIdxAbove = find(sum(rPDFF>=InflP,2));
    KeptIdxBelow = find(sum(rPDFF<InflP,2));
end

NoiseCorrResults.r_plus_cells = KeptIdxBelow;
NoiseCorrResults.r_minus_cells = KeptIdxAbove;

%%%%%%%% Plot change in noise correlation for < and >=InflP %%%%%%%%
set(0,'CurrentFigure',Fig) 

subplot(2,2,1)
y_Lims = get(gca,'Ylim')
plot([InflP InflP] ,[y_Lims(1) y_Lims(2)],'k--' )

subplot(2,2,2)
%>=InflP
rdiff_minus=rHDFF(KeptIdxAbove) - rPDFF(KeptIdxAbove);
mu = nanmean(rdiff_minus);
%bins = linspace(min(rdiff),max(rdiff),20);
CorrHist = hist(rdiff_minus,bins);
CorrHist = CorrHist./length(rPDFF);
bar(bins,CorrHist,1,'facecolor',[.33 .33 .33],'edgecolor','w')
NoiseCorrResults.r_minus = rdiff_minus;

try
aa=axis;
hold on
plot([mu mu],[0 aa(4)],'k--')
%<InflP
rdiff_plus=rHDFF(KeptIdxBelow) - rPDFF(KeptIdxBelow);
NoiseCorrResults.r_plus = rdiff_plus;
mu = nanmean(rdiff_plus);
%bins = linspace(min(rdiff),max(rdiff),20);
CorrHist = hist(rdiff_plus,bins);
CorrHist = CorrHist./length(rPDFF);
bar(bins,CorrHist,1,'facecolor',[.67 .67 .67],'edgecolor','w')
aa=axis;
hold on
plot([mu mu],[0 aa(4)],'k--')
ylim([0 .15])
xlim([- .6 .6])
plot([0 0],[0 aa(4)],'k')
xlabel('\DeltaNoise Corr.')
ylabel('Probability')
set(gca,'fontsize',8)
title('\DeltaNoise Correlation (r)')
catch
end 
%%%%%%%% Statisitical tests for hit-pass about InflP %%%%%%%% 
% %>=InflP
% A = rHDFF(KeptIdxAbove)-rPDFF(KeptIdxAbove);
% [h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
% NoiseCorrResults.HitPassAbove.mu=nanmean(A);
% NoiseCorrResults.HitPassAbove.std=2*std(A)./sqrt(length(A));
% NoiseCorrResults.HitPassAbove.bootsig001=[h p];
% %<InflP
% KeptIdx = find(sum(rPDFF<InflP,2));
% A = rHDFF(KeptIdxBelow)-rPDFF(KeptIdxBelow);
% [h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
% NoiseCorrResults.HitPassBelow.mu=nanmean(A);
% NoiseCorrResults.HitPassBelow.std=2*std(A)./sqrt(length(A));
% NoiseCorrResults.HitPassBelow.bootsig001=[h p];
% %< vs. >= InflP
% A=[];
% A = abs(rHDFF(KeptIdxAbove)-rPDFF(KeptIdxAbove));
% B=[];
% B = abs(rHDFF(KeptIdxBelow)-rPDFF(KeptIdxBelow));
% [h p]=bootsig(A,B, 10000, 0.001)
% NoiseCorrResults.HitPassBelowAbove.mu=nanmean(A)-nanmean(B);
% [h p that]=bootsig(A,B, 10000, 0.001);
% NoiseCorrResults.HitPassBelowAbove.std=2*std(that);
% NoiseCorrResults.HitPassBelowAbove.bootsig=[h p];


% Percentage r+ / r- Analysis 
% 
% n_neurons = length(rHDFF);
% bar( [length(KeptIdxAbove), length(KeptIdxBelow)]/ n_neurons * 100 )
% ylabel('%')
% xticks(1:2)
% xticklabels({'r+','r-'})
% xlabel('\Deltar')

%Analysis: Active-Passive for Tar/BF Octave diff for cell pairs with >=InflP correlations
F1idx_above = intersect(find(close_idx),KeptIdxAbove);
F2idx_above = intersect(find(~close_idx),KeptIdxAbove);
F1rPDFF_a = rPDFF(F1idx_above);
F1rHDFF_a= rHDFF(F1idx_above);
F2rPDFF_a= rPDFF(F2idx_above);
F2rHDFF_a= rHDFF(F2idx_above);

F1idx_below = intersect(find(close_idx),KeptIdxBelow);
F2idx_below = intersect(find(~close_idx),KeptIdxBelow);
F1rPDFF_b = rPDFF(F1idx_below);
F1rHDFF_b= rHDFF(F1idx_below);
F2rPDFF_b= rPDFF(F2idx_below);
F2rHDFF_b= rHDFF(F2idx_below);

subplot(2,2,3) 
NikTarBFPlot(F1rPDFF_a, F1rHDFF_a,F2rPDFF_a, F2rHDFF_a,PlotColor)
title('\DeltaTar (octaves)')
subplot(2,2,4) 
NikTarBFPlot(F1rPDFF_b, F1rHDFF_b,F2rPDFF_b, F2rHDFF_b,PlotColor)


% Above inF = r- 
% below inf = r+
NoiseCorrResults.BF.Close.Passive.Above = F1rPDFF_a;
NoiseCorrResults.BF.Close.Active.Above = F1rHDFF_a;
NoiseCorrResults.BF.Close.Passive.Below = F1rPDFF_b;
NoiseCorrResults.BF.Close.Active.Below = F1rHDFF_b;

NoiseCorrResults.BF.Far.Passive.Above = F2rPDFF_a;
NoiseCorrResults.BF.Far.Active.Above = F2rHDFF_a;
NoiseCorrResults.BF.Far.Passive.Below = F2rPDFF_b;
NoiseCorrResults.BF.Far.Active.Below = F2rHDFF_b;





%Plot active & passive


function NikTarBFPlot(F1rPDFF, F1rHDFF,F2rPDFF, F2rHDFF,PlotColor)
    
%Plot active & passive
mu = [nanmean(F1rPDFF); nanmean(F1rHDFF)];
ste = [2.*nanstd(F1rPDFF)/sqrt(length(F1rPDFF)); 2.*nanstd(F1rHDFF)/sqrt(length(F1rHDFF))];
errorbar(1:2,mu,ste,PlotColor,'linewidth',2)
hold on
mu = [nanmean(F2rPDFF); nanmean(F2rHDFF)];
ste = [2.*nanstd(F2rPDFF)/sqrt(length(F2rPDFF)); 2.*nanstd(F2rHDFF)/sqrt(length(F2rHDFF))];
errorbar(1.3:2.3,mu,ste, [PlotColor '--'],'linewidth',2)
axis tight
xlim([0 3])
aa=axis;
%ylim([aa(3)-.05 aa(4)+.05])
ylim([ -.15 .45])
plot([aa(1) aa(2)],[0 0],'k')
ylabel('r')
set(gca,'xtick',1:2)
set(gca,'xticklabel',{'Passive','Hit'})
set(gca,'fontsize',10)

 
function  [h p that]= bootsig(x,y, nexps, alpha, N, method, ciplot)

if nargin < 3
    nexps = 10000;
    alpha = 0.05;
    N=[];
    method = 'count';
    ciplot = 0;
end
if  nargin >= 3 && nargin < 4
    alpha = 0.05;
    N=[];
    method = 'count';
    ciplot = 0;
end
if  nargin >= 4 && nargin < 5
    N=[];
    method = 'count';
    ciplot = 0;
end
if  nargin >= 5 && nargin < 6
    method = 'count';
    ciplot = 0;
end
if  nargin >= 6 && nargin < 7
    ciplot = 0;
end
if isempty(y)
    y=zeros(size(x));
end
if iscell(x) && iscell(y)
    if isempty(N)
        L = min([length(x{1}) length(x{2}) length(y{1}) length(y{2})]);
    else
        L = N;
    end
else
    if isempty(N)
        L = min([length(x) length(y)]);
    else
        L = N;
    end
end
if strcmpi(method,'count')
    t = abs(nanmean(y) - nanmean(x));
    null = [x; y];
    that=[];
    for i =1:nexps
        xhat = randsample(null,L,'true');
        yhat = randsample(null,L,'true');
        that = [that; abs(nanmean(yhat) - nanmean(xhat))];
    end
    h=0;
    p = sum(that>t)./nexps;
    if p < alpha
        h = 1;
    end
elseif strcmpi(method,'ci')
    if iscell(x) && iscell(y)
        that=[];
        for i =1:nexps
            x1 = randsample(x{1},L,'true');
            x2 = randsample(x{2},L,'true');
            xdiff = x1-x2;
            y1 = randsample(y{1},L,'true');
            y2 = randsample(y{2},L,'true');
            ydiff=y1-y2;
            that = [that; abs(nanmean(ydiff) - nanmean(xdiff))];
        end
        ste = std(that);
    else
        t = (nanmean(y) - nanmean(x));
        null = [x; y];
        
        that=[];
        for i =1:nexps
            xhat = randsample(null,L,'true');
            yhat = randsample(null,L,'true');
            that = [that; (nanmean(yhat) - nanmean(xhat))];
        end
    end
    p = alpha;
    h=0;
    ci = [prctile(that,100*alpha/2) prctile(that,100-(100*alpha/2))];
    if ci(1,1)*ci(1,2) > 0
        h=1;
        
    end
    if ciplot
        figure
        ci = [prctile(that,100*alpha/2) prctile(that,100-(100*alpha/2))];
        hist(that,50)
        hold on
        plot([ci(1,1) ci(1,1)],[0 700],'g','linewidth',3)
        plot([ci(1,2) ci(1,2)],[0 700],'g','linewidth',3)
        plot([t t],[0 700],'r')
        
    end
end

function  [h p]= bootsig0(x, nexps, alpha, ciplot)

if nargin < 3
    alpha = 0.05;
    ciplot = 0;
    
elseif nargin < 4
    ciplot = 0;
    
end
p=alpha;
that=[];
parfor i =1:nexps
    xhat = randsample(x,length(x),'true');
    that = [that; mean(xhat)];
    
end

h=0;
t=0; %Test statstic, ie. is the mean sig. diff than t=0?
ci = [prctile(that,100*alpha/2) prctile(that,100-(100*alpha/2))];
if mean(x) < 0
    if t  > ci(2)
        h=1;
    end
    
elseif mean(x) > 0
    if t <ci(1)
        h=1;
    end
    
end

if ciplot
    figure
    ci = [prctile(that,100*alpha/2) prctile(that,100-(100*alpha/2))];
    hist(that,50)
    hold on
    plot([ci(1,1) ci(1,1)],[0 700],'g','linewidth',3)
    plot([ci(1,2) ci(1,2)],[0 700],'g','linewidth',3)
    plot([t t],[0 700],'r')
end



    function [out] = getSubtractedStats(Passive,Active)
        
        
        
        out = struct();
        gain = Active - Passive;
        mean = nanmean(gain);
        
        var1 = nanvar(Passive);
        n1 = length(Passive);
        var2 = nanvar(Active);
        n2 = length(Active);
        
        pooled_var = (var1 * (n1 -1) + var2 * (n2-1))  / (n1 + n2 - 2);
        
        CI = sqrt(pooled_var / (n1+n2) ) * 1.96;
        
        out.passive = Passive;
        out.active = Active;
        out.gain = gain;
        out.mean = mean;
        out.var = pooled_var;
        out.n = n1+n2;
        out.CI = CI;