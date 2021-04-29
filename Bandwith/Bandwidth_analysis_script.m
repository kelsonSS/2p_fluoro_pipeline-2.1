yp = BandwidthAnalysis(FRA_young_noise, 'BRFS',0,'Pos',.5);
ap = BandwidthAnalysis(FRA_old_noise, 'BRFS',0,'Pos',.5);

an =  BandwidthAnalysis(FRA_old_noise, 'BRFS',0,'Neg',.5);
yn = BandwidthAnalysis(FRA_young_noise, 'BRFS',0,'Neg',.5);

an = an{1};
ap = ap{1};
yn = yn{1};
yp = yp{1};


all_p = cat(2,ap,yp);
all_n = cat(2,an,yn);

all_p = all_p(:);
all_n = all_n(:);



old_idx = repmat({'old'}, numel(ap),1);
young_idx = repmat({'young'},numel(yp),1);
age_idx = cat(1,old_idx,young_idx);

clear yn yp an ap

lvl_idx = repmat({'99';'20';'10';'0'},1,length(all_p)/4) ;
lvl_idx = lvl_idx(:);

[x,y,stats,z]= anovan(all_p,{age_idx,lvl_idx},'model','interaction');
ap_stats = multcompare(stats,'Dimension',[1 2]);
ap_stats = ap_stats(ap_stats(:,1)+1 ==ap_stats(:,2),:);
ap_stats = ap_stats(1:2:end,:);


[x,y,stats,z]= anovan(all_n,{age_idx,lvl_idx},'model','interaction');
an_stats= multcompare(stats,'Dimension',[1 2]);
an_stats = an_stats(an_stats(:,1)+1 ==an_stats(:,2),:);
an_stats = an_stats(1:2:end,:)

clear x y z age_idx young_idx old_idx all_n all_p



