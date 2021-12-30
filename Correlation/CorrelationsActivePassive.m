function out = CorrelationsActivePassive(TNBehavior,SaveName,BehaviorType,PassiveFreqs)

if ~exist('SaveName','var')
    SaveName = [];
end 

if ~exist('BehaviorType','var')
    BehaviorType = 'Hit'   % What subset of behavior to select: Hit,Miss,Early
end 



if ~exist('PassiveFreqs','var')
    PassiveFreqs = 11314;
end 

expts = size(TNBehavior,1);
passive_all = [];
passive_all2 =[];
active_all2 =[];
active_all = [];
corrs_passive_expt = {};
corrs_active_expt = {};
corrs_passive_expt2 = {};
corrs_active_expt2 = {};

test_passive = CorrelationsByAnimalBehavior(TNBehavior(:,1));
test_active = CorrelationsByAnimalBehavior(TNBehavior(:,2));

[passive_NCorr,~] = MungeCorr2Mat(test_passive.NCorr);
[active_NCorr, ~] = MungeCorr2Mat(test_active.NCorr);

figure;
hold on
cdfplot(passive_NCorr)
cdfplot(active_NCorr)
title('Noise Correlations Old Method CDF')
if SaveName
saveas(gcf, sprintf('%s_NoiseCorrelations_K.pdf',SaveName))
end 

legend({'Passive','Active'})

for expt = 1:expts
    
    passive = TNBehavior{expt,1} ;
    active =  TNBehavior{expt,2};
    
   
    excited = squeeze(max(nanmean(active.DFF(30:end,:,:)))) > 0;
    %passive_idx = passive.active{:,2}>0; 
    %active_idx  = active.active{:,2}>0;
    passive_idx = passive.responsive;
    active_idx = active.responsive;
    
    uLevels_expt =sort(unique(active.FreqLevelOrder{:,2}),'descend');
    
    if length(passive_idx) ~= length(active_idx) || sum( active_idx & passive_idx) == 0
        sprintf('%d ',expt) 
        continue
    end 
    responsive_idx = squeeze(passive_idx & active_idx & excited)
    responsive = table([1:length(passive_idx)]',responsive_idx,'VariableNames',{'Neuron','Activity'});
    
    passive.active = responsive;
    active.active = responsive;
   
    switch BehaviorType
        case 'Hit'
            behave_idx = find(active.handles{1}.Hits);
          
        case 'Miss'
              behave_idx = find(active.handles{1}.Miss);
        
        case 'Early'
            behave_idx = find(~active.handles{1}.Early);
              
        case 'NotHit'
            behave_idx = find(~active.handles{1}.Hits);
        otherwise 
            error('Unrecognized behave index')
    end
    behave_idx = behave_idx( behave_idx<= size(active.DFF,2));
    
    active.DFF = active.DFF(:,behave_idx,:);
    active.DFF_Z = active.DFF_Z(:,behave_idx,:);
    active.FreqLevelOrder = active.FreqLevelOrder(behave_idx,:);
        
    
    passive_freq_idx = any(passive.FreqLevelOrder{:,1} == PassiveFreqs,2);
    passive.DFF = passive.DFF(:,passive_freq_idx,:);
    passive.DFF_Z = passive.DFF_Z(:,passive_freq_idx,:);
    passive.FreqLevelOrder = passive.FreqLevelOrder(passive_freq_idx,:);
    passive.BestFrequencies = BestFrequencyAnalysis(passive);

    
     TNBehavior{expt,1} = passive ;
     TNBehavior{expt,2} = active;
end 

test_passive2 = CorrelationsByAnimalBehavior(TNBehavior(:,1));
test_active2 = CorrelationsByAnimalBehavior(TNBehavior(:,2));
[passive_all2,passive_all2_mu] = MungeCorr2Mat(test_passive2.NCorr);
[active_all2, active_all2_mu] = MungeCorr2Mat(test_active2.NCorr);


  %  corrs_passive = CorrelationsByAnimal(passive);
  %  corrs_active = CorrelationsByAnimal(active);
    %catch
    %    continue
    %end
   
%     corrs_passive_expt{expt} = corrs_passive.NCorr{1};
%     corrs_passive_expt2{expt} = corrs_passive.NCorrCell{1};
%     corrs_active_expt{expt} = corrs_active.NCorr;
%     corrs_active_expt2{expt} = corrs_active.NCorrCell;
%     
%     
%     passive_all = cat(1,passive_all,corrs_passive.NCorr{1}(:));
%     passive_all2 = cat(1,passive_all,corrs_passive.NCorrCell{1}(:));
%     
%     active_corrs_levels = subsetActiveCorrs(corrs_passive.NCorr,corrs_active.NCorr,SNRs,uLevels_expt);
%     active_all = cat(1,active_all,active_corrs_levels);
%     
%     active_all2 = cat(1,active_all2,corrs_active.NCorrCell{1}(:));
%     
       
%  end
if SaveName
CorrsByDistance(TNBehavior(:,1),10,[SaveName '_Passive'])
CorrsByDistance(TNBehavior(:,2),10,[SaveName '_Active'])
end 


close_idx_all  = CompareBestFrequencyToBehaviorTarget(TNBehavior);

close_idx = cell2mat(close_idx_all);

out = NikCompareCorrs(passive_all2_mu,active_all2_mu,close_idx);
% figure;
% hold on
% cdfplot(passive_all2)
% cdfplot(active_all2)
% legend({'Passive','Active'})
% title('Noise Correlations')
if SaveName
saveas(gcf, sprintf('%s_NoiseCorrelations_Nik.pdf',SaveName))
end 







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



    function  NoiseCorrResults = NikCompareCorrs(rPDFF, rHDFF,close_idx)

NoiseCorrResults = [];

Fig = figure

subplot(2,2,1)

n = min(length(rHDFF),length(rPDFF));

rPDFF = rPDFF(1:n);
rHDFF = rHDFF(1:n);


bins = linspace(min(rPDFF),max(rPDFF),20);
%Passive
CorrHist = hist(rPDFF,bins);
CorrHist = CorrHist./sum(CorrHist);
plot(bins,cumsum(CorrHist),'color',[.33 .33 .33],'linewidth',2)
hold on
%Hit
bins = linspace(min(rHDFF),max(rHDFF),20);
CorrHist = hist(rHDFF,bins);
CorrHist = CorrHist./sum(CorrHist);
plot(bins,cumsum(CorrHist),'b','linewidth',2)
xlim([-.6 .6])
plot([0 0],[0 1],'k')
xlabel('Noise Corr.')
ylabel('Probability')
set(gca,'fontsize',8)
title('Noise Correlation (r) CDF')

%%%%%%%% Statisitical tests for hit vs pass %%%%%%%% 
%hit
A = rHDFF;
[h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
NoiseCorrResults.Hit.mu=nanmean(A(~isnan(A)));
NoiseCorrResults.Hit.std=2*std(A(~isnan(A)))./sqrt(length((~isnan(A))));
NoiseCorrResults.Hit.bootsig001=[h p];
%passive
A = rPDFF;
[h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
NoiseCorrResults.Pass.mu=nanmean(A(~isnan(A)));
NoiseCorrResults.Pass.std=2*std(A(~isnan(A)))./sqrt(length((~isnan(A))));
NoiseCorrResults.Pass.bootsig001=[h p];
%Hit vs passive
A = rHDFF;
B = rPDFF;
[h p]=bootsig(A(~isnan(A)), B(~isnan(B)), 10000, 0.001);
NoiseCorrResults.HitPass.mu=nanmean(A(~isnan(A)))-nanmean(B(~isnan(B)));
NoiseCorrResults.HitPass.std=2*nanstd(A-B)./sqrt(length(A-B));
NoiseCorrResults.HitPass.bootsig=[h p];

%%%%%%%% Find Hit-Passive CDF diff inflection point
CorrHist = hist(rPDFF,bins);
PCorrHist = cumsum(CorrHist./sum(CorrHist));
CorrHist = hist(rHDFF,bins);
HCorrHist = cumsum(CorrHist./sum(CorrHist));
diffCorrHist = HCorrHist-PCorrHist;
[junk m]=max(diff(diffCorrHist));
InflP = bins(m);
NoiseCorrResults.InflP = InflP;
KeptIdxAbove = find(sum(rPDFF>=InflP,2));
KeptIdxBelow = find(sum(rPDFF<InflP,2));

%%%%%%%% Plot change in noise correlation for < and >=InflP %%%%%%%%
set(0,'CurrentFigure',Fig) 
subplot(2,2,2)
%>=InflP
rdiff=rHDFF(KeptIdxAbove) - rPDFF(KeptIdxAbove);
mu = nanmean(rdiff);
bins = linspace(min(rdiff),max(rdiff),20);
CorrHist = hist(rdiff,bins);
CorrHist = CorrHist./length(rPDFF);
bar(bins,CorrHist,1,'facecolor',[.33 .33 .33],'edgecolor','w')

aa=axis;
hold on
plot([mu mu],[0 aa(4)],'k--')
%<InflP
rdiff=rHDFF(KeptIdxBelow) - rPDFF(KeptIdxBelow);
mu = nanmean(rdiff);
bins = linspace(min(rdiff),max(rdiff),20);
CorrHist = hist(rdiff,bins);
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

%%%%%%%% Statisitical tests for hit-pass about InflP %%%%%%%% 
%>=InflP
A = rHDFF(KeptIdxAbove)-rPDFF(KeptIdxAbove);
[h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
NoiseCorrResults.HitPassAbove.mu=nanmean(A);
NoiseCorrResults.HitPassAbove.std=2*std(A)./sqrt(length(A));
NoiseCorrResults.HitPassAbove.bootsig001=[h p];
%<InflP
KeptIdx = find(sum(rPDFF<InflP,2));
A = rHDFF(KeptIdxBelow)-rPDFF(KeptIdxBelow);
[h p]=bootsig0(A(~isnan(A)), 10000, 0.001);
NoiseCorrResults.HitPassBelow.mu=nanmean(A);
NoiseCorrResults.HitPassBelow.std=2*std(A)./sqrt(length(A));
NoiseCorrResults.HitPassBelow.bootsig001=[h p];
%< vs. >= InflP
A=[];
A = abs(rHDFF(KeptIdxAbove)-rPDFF(KeptIdxAbove));
B=[];
B = abs(rHDFF(KeptIdxBelow)-rPDFF(KeptIdxBelow));
[h p]=bootsig(A,B, 10000, 0.001)
NoiseCorrResults.HitPassBelowAbove.mu=nanmean(A)-nanmean(B);
[h p that]=bootsig(A,B, 10000, 0.001);
NoiseCorrResults.HitPassBelowAbove.std=2*std(that);
NoiseCorrResults.HitPassBelowAbove.bootsig=[h p];


% Percentage r+ / r- Analysis 
subplot(2,2,3) 
n_neurons = length(rHDFF);
bar( [length(KeptIdxAbove), length(KeptIdxBelow)]/ n_neurons * 100 )
ylabel('%')
xticks(1:2)
xticklabels({'r+','r-'})
xlabel('\Deltar')

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




subplot(2,2,4)
hold on 
NikTarBFPlot(F1rPDFF_a, F1rHDFF_a,F2rPDFF_a, F2rHDFF_a)
NikTarBFPlot(F1rPDFF_b, F1rHDFF_b,F2rPDFF_b, F2rHDFF_b)
ylim([ -.15 .35])

%Plot active & passive


function NikTarBFPlot(F1rPDFF, F1rHDFF,F2rPDFF, F2rHDFF)
    
%Plot active & passive
mu = [nanmean(F1rPDFF); nanmean(F1rHDFF)];
ste = [2.*nanstd(F1rPDFF)/sqrt(length(F1rPDFF)); 2.*nanstd(F1rHDFF)/sqrt(length(F1rHDFF))];
errorbar(1:2,mu,ste,'k','linewidth',2)
hold on
mu = [nanmean(F2rPDFF); nanmean(F2rHDFF)];
ste = [2.*nanstd(F2rPDFF)/sqrt(length(F2rPDFF)); 2.*nanstd(F2rHDFF)/sqrt(length(F2rHDFF))];
errorbar(1:2,mu,ste,'k--','linewidth',2)
axis tight
xlim([0 3])
aa=axis;
ylim([aa(3)-.05 aa(4)+.05])
plot([aa(1) aa(2)],[0 0],'k')
ylabel('r')
set(gca,'xtick',1:2)
set(gca,'xticklabel',{'Passive','Hit'})
set(gca,'fontsize',10)
title('\DeltaTar (octaves)')
 
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



