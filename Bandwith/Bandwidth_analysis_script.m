
% entropy analysis 
H = cellfun(@(x) x.entropy,AVGBF,'UniformOutput',0);

for ii = 1:3
HH{ii} = cell2mat(H(:,ii));
end 


H_mu = cellfun(@mean, HH);
H_std = cellfun(@std,HH) ;



figure
hold on 
bar(H_mu)
errorbar([1,2,3],H_mu,H_std/sqrt(6),'.')

% Anova to determine if Entropy is different among neurons  
[p,~,sts]=anova1(H,[],'off');
H_stats = multcompare(sts);


%% 2D-Bandwith Analysis



DFL = cellfun(@(x) x.Lvls, AVGBF,'UniformOutput',0); 

% bandwith grouping
% using bandwith calculation but using entire 2D dataset. (Currently: No
% Tuning assumed)
BD = cellfun(@Band2D, DFL,'UniformOutput',0);


BD_class =  { cell2mat(BD(:,1)'),...
              cell2mat(BD(:,2)'),...
              cell2mat(BD(:,3)') };
                   
          
BD_mu = cellfun(@mean, BD_class);


BD_ste = cellfun(@(x) std(x)/sqrt(length(x)), BD_class);

% bandwith by level
BDL = cellfun(@(x) Band2D(x,4) , DFL,'UniformOutput',0);

BDL_class =  { cell2mat(BDL(:,1)'),...
              cell2mat(BDL(:,2)'),...
              cell2mat(BDL(:,3)') };
                   
          
BDL_mu = cell2mat(cellfun(@(x) mean(x,2),...
                  BDL_class,'UniformOutput',0));


BDL_ste = cell2mat(cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),...
                   BDL_class,'UniformOutput',0));





%plotting Bandwidth
figure
axis square
hold on 
bar(BD_mu,'k')
errorbar([1,2,3],BD_mu,BD_ste,'k.','LineWidth',3)
xticks(1:3)
xticklabels({'Tones','Noise','Offset'})
xlabel('Class')
ylabel(' 2D-Bandwidth (AU)')


% plotting Bandwith by level
figure 
axis square
errorbar(BDL_mu,BDL_ste,'LineWidth', 2)
set(gca,'Xdir','reverse')
xlim([0 5])
xticks([0:5])
xticklabels({[],'+10','+20','+30','inf'})
xlabel('SNR (DB)')
ylabel('Bandwith (AU)')
legend('Tones','Noise','Offset')



%% Signficance testing 
% grouping into list of all BDs and list of group IDs
class = [];
BD_mu_all= []; 
for ii = 1:length(BD_class)
    BD_mu_all = vertcat(BD_mu_all, BD_class{ii}'); 
   
    class = vertcat(class, ones( length(BD_class{ii}) , 1  ) * ii );
    
end
    
% Signficance testing 
[BD_p,~,BD_sts] = anovan(BD_mu_all,class);
[BDS.comp,BDS.sig,BDS.H,BDS.nms] = multcompare(BD_sts);    



% 2D signficance testing 
BDL_expt = cell(size(BDL));
BDL_lvl =  cell(size(BDL));
BDL_class_ID =  cell(size(BDL));
for expt = 1:6
    for grp = 1:3
        [r,c] = size(BDL{expt,grp}); 
       
    BDL_expt{expt,grp} = ones(r,c)* expt ;
    BDL_lvl{expt,grp}  = repmat([1:4]',1,c);
    BDL_class_ID{expt,grp} = ones(r,c) * grp;
    
    end 
end 


% Bandwith by Level

[BDL_p,~,BDL_sts]= anovan(CellFlat2D(BDL),{...
                        CellFlat2D(BDL_lvl),...
                        CellFlat2D(BDL_Class_ID)},'model','interaction',...
                       'varnames',{'Level','Class'});




[BDLS.comp,BDLS.sig,BDLS.h,BDLS.nms] = multcompare(BDL_sts,'dimension',[1,2]);



BDL_Mu = cellfun(@(x) mean(x,2), BDL,'UniformOutput',0));
















function flat = CellFlat2D(C)

mat = [C{:,:}];
flat = mat(:);

end 