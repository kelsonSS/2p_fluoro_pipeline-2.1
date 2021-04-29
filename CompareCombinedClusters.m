function p = CompareCombinedClusters(young,old)


old_clust = MungeClusters(old);
young_clust = MungeClusters(young);

%
[~,p] = ttest2(old_clust.Diversity,young_clust.Diversity);

%plotting
figure
errorbar(2,mean(old_clust.Diversity),std(old_clust.Diversity) / sqrt(12))
errorbar(1,mean(young_clust.Diversity),std(young_clust.Diversity) / sqrt(12))
hold on
errorbar(2,mean(old_clust.Diversity),std(old_clust.Diversity) / sqrt(12))
xlim([0,4])
xlabel('Cluster Diversity')
ylim([0,6])
legend({'young','old'})



end 


function Out = MungeClusters(c)

clust = c.Combined_Classes;
clust_idx = clust ~= 0;
expt_id = c.experiment_list(clust_idx);

% munge clusters
% remove inactive neurons 
clust(clust == 0) = [];
% sort clusters 
clust(clust ==1) = 10;
clust = clust - 1; 

% merge clusters with expt ids
munged = cat(2,clust,expt_id);


n = max(munged(:,1)); % n clusters
m = max(munged(:,2)); % m experiments 


Out = [];
Out.Clust = munged;
Out.Population_counts = histcounts(munged(:,1),n);
Out.All =[];
Out.Diversity = [];
Out.Pos = [];
Out.Neg = [];
Out.Mixed = [];

%WARNING- these values are hard coded for this expt!! only use ALL vals if not
%using for passive expt!!!
out_id = 1;
for expt = 1:m
   
    expt_idx = munged(:,2) == expt;
    
    %if sum(expt_idx) < 20
    %    continue
    % end 
    
    s = munged(expt_idx,1);
    Out.All(:,out_id) = histcounts(s,n) ./ length(s);
    
    Out.Pos(out_id) = sum( s <=4 ) ./length(s) ;
    Out.Neg(out_id) = sum( s<=7 & s >4) ./ length(s);
    Out.Mixed(out_id) = sum( s>7)./length(s);
    
    out_id = out_id+1;
  
    
    
    
    
    

    
    
end 

Out.Diversity = sum(Out.All >.05);



end 

function chi = ChiSquareTest(Population)

n = sum(Population);
m = length(Population);
expected  = n/m;

chi = sum( [Population - expected].^2/expected );

end 

function p = Bootstrap(sample,k,target_val)

% k = number of repeats
s = 200; % sample_size

x =  reshape(randsample(sample,k*s,true),k,1000);

x = sum(x,2)./ s;

p =  sum( x > target_val) / length(x);


end 


function out = Bootstrap_all(young_clust,old_clust)


old_neg = mean(old_clust.Neg);
old_pos = mean(old_clust.Pos);
old_mixed = mean(old_clust.Mixed);

yc = young_clust.Clust(:,1);
ac = old_clust.Clust(:,1);

out = []
out.Neg_P = Bootstrap(createBootstrapSample(yc,ac,'Neg') ,...
                  10000,old_neg);
out.Pos_P = Bootstrap(createBootstrapSample(yc,ac,'Pos') ,...
                  10000,old_pos);
out.Mixed_P = Bootstrap(createBootstrapSample(yc,ac,'Mixed') ,...
                  10000,old_mixed);
              
              
              
end 
       
function out = createBootstrapSample(sample_1,sample_2,type)
% creates a mixed sample as a null distribution


s = cat(1,sample_1,sample_2)

switch type
    case 'Neg'
   out = s >4 & s <=7; 
    case 'Pos'
        out = s <=4;
    case 'Mixed'
        out = s<7; 
end 
    
end 
    
    










