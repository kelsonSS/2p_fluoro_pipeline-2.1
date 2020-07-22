function PlotAvgTracePassive(DFF)

% takes a time X trial X neuron normalized matrix or a  
 %        time X Neuron  matrix created by averaging the first matrix 
 % over trials and outputs the avgtrace figure as seen in the passivepaper. 
 % currently support for 90frame trials and 150 frame trials. assumes 30fps 
if ndims(DFF)<2 

elseif ndims(DFF) == 2
    n_neurons = size(DFF,2);
elseif ndims(DFF) == 3
    n_neurons = size(DFF,3);
    DFF = reshape(DFF,size(DFF,1) ,[]);
end 

mu = squeeze(nanmean(DFF,2));
CI = squeeze(nanstd(DFF,[],2)) / sqrt(n_neurons) * 1.96 ;
%    


% plotting 
figure
shadedErrorBar([],mu,CI);

title( sprintf('%d Neurons', n_neurons ));
ylim([-1 1])
aa = axis;
hold on

plot([aa(1) aa(2)], [0 0 ] ,'k--')
if size(DFF,1) >90
    plot([30 30], [aa(3) aa(4)],'r--')
    plot([60 60], [aa(3) aa(4)],'g--')
    plot([90 90], [aa(3) aa(4)],'g--')
    plot([120 120], [aa(3) aa(4)],'r--')
    axis tight
    xlim([0 150])
else
    plot([30 30], [aa(3) aa(4)],'g--')
    plot([60 60], [aa(3) aa(4)],'g--')
    axis tight
    xlim([0 90])
end
        
