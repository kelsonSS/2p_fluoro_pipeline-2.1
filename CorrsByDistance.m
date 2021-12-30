function Out = CorrsByDistance(DF,BinSizeMicrons,SaveName)

if ~exist('SaveName','var')
    SaveName = ''
end 

try
corrs_all = DF.CorrByAnimal;
catch
corrs_all = CorrelationsByAnimalBehavior(DF);
end 

[CellDists,max_dist]= getCellDistances(DF);


% empirically there are too few cells that have >400uM distance on B-scope to
% have meaningful comparisions so we cut off there
max_dist = min(max_dist,300);

C= struct();
fn = fieldnames(corrs_all);
corr_idx = ( contains(fn,'LCorr') | contains(fn,'NCorr') ) &...
        ( ~( contains(fn,'Animal')| contains(fn,'Diff')|  contains(fn,'Cell') )) ;
fn = fn(corr_idx);



Out.BinSize = BinSizeMicrons;
Out.max_dist_condisered = max_dist;

   for field_idx = 1:length(fn)
       try
       Out.(fn{field_idx}) = getCorrsAtDistance(corrs_all.(fn{field_idx}),...
                                           CellDists,max_dist,BinSizeMicrons);

       PlotCorrsDistance(Out.(fn{field_idx}), [SaveName '-' fn{field_idx}], BinSizeMicrons)
       catch
       end 

   end 

 







function CorrDist = getCorrsAtDistance(Corrs,Distances,MaxDist,BinSize)

CorrDist.Cell = {};
CorrDist.Animal = {};

insertion_idx = 1;
    for DistanceLowerBound = 1:BinSize:MaxDist

        DistanceUpperBound = DistanceLowerBound+BinSize;
         
      [CorrDist.Cell{insertion_idx},CorrDist.Cell_mu(insertion_idx,:),CorrDist.Cell_CI(insertion_idx,:),...
       CorrDist.Animal_mu(insertion_idx,:),CorrDist.Animal_CI(insertion_idx,:) ] = ...
                        subsetCorrsbyDistance(Corrs,Distances,...
                                              DistanceLowerBound,...
                                              DistanceUpperBound);
    
    
        insertion_idx = insertion_idx+1;
    
    end 



function [Cell,Cell_mu,Cell_CI,Animal_mu,Animal_CI] = subsetCorrsbyDistance(C,dists,lower,upper)

Cell = [];
Animal = {};

% extract
for expt = 1: length(dists)

      corrs_animal = C(expt,:);
      dist_idx = find( dists{expt} >= lower & dists{expt} < upper)  ;
      corrs_animal = cellfun(@(x) x(dist_idx),corrs_animal,'UniformOutput',0);
        
      Animal(expt,:) = corrs_animal;
      Cell = cat(1,Cell, cell2mat(corrs_animal));
  




end 

% process
Cell_mu = nanmean(Cell); 
Cell_CI = nanstd(Cell) / sqrt( length(Cell)) * 3.291; %99.9 CI

Animal = cellfun(@nanmean,Animal);
Animal_mu = nanmean(Animal); 
Animal_CI = nanstd(Animal) / sqrt( length(Animal)) * 1.96;

 

function  PlotCorrsDistance(CellCorrs,FieldName, BinSizeMicrons)


x = [0:BinSizeMicrons:(size(CellCorrs.Cell_mu,1)-1)*BinSizeMicrons];
figure
set(gcf,'Position',[660 450 560 520])


if size(CellCorrs.Cell_mu,2) == 1
    
    shadedErrorBar(x,CellCorrs.Cell_mu,CellCorrs.Cell_CI)
    title('All')
    stylefig()
else
    levelTitles = {'Tones','+20','+10','0'}
    
    for ii =1:4
        subplot(2,2,ii)
        shadedErrorBar(x,CellCorrs.Cell_mu(:,ii),CellCorrs.Cell_CI(:,ii))
        title(levelTitles{ii})
        stylefig()
        
    end
end
 
suptitle(FieldName)

saveas(gcf,sprintf('%s-ByDistance.pdf', FieldName))


function stylefig()
    ylim([0 .5])
    ylabel('Correlation (AU)')
    xlabel('Distance (um)')







        

    
 






