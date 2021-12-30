function plotAllDFFs(DFF,norm,classes)
    
if ~exist('classes','var')
         classes = false; 
         num_classes = 1;
        
end 

if ~exist('norm','var')
    norm = 1;
end 
    

if isnumeric(DFF) && ~isstruct(DFF)
    
    plotDFFs(DFF,'All',size(DFF,3))
else 
    if classes 
         num_classes = max(DFF.Combined_Classes);
         %ClassNames = DFF.Classes;
         
         % get classes of all active neurons
     All_Classes =  DFF.Combined_Classes .* (DFF.active{:,2}>0);  
         
         
     for ii = 1:num_classes 
         Class_idx{ii} =  All_Classes == ii ;
     end 
    else
        num_classes = 1;
        Class_idx{1} = DFF.active{:,2}>0; 
        ClassNames{1} = 'All';
    end
         
         
    for class = 1:num_classes
        
        %% Selection
        
        if isfield(DFF,'Clean_idx')
            Clean_idx = DFF.Clean_idx;
        elseif isfield(DFF,'ArtifactIndex')
            Clean_idx = ~DFF.ArtifactIndex;
        end
        
        final_idx = Clean_idx & Class_idx{class}
           
        
        if norm
            try 
            Fluoro = DFF.DFF_Z(:,:,final_idx);
            % baseline correction
            B_Vec = repmat(nanmean(Fluoro(1:30,:,:)),[size(Fluoro,1),1,1]);
            Fluoro = Fluoro - B_Vec;
            
            catch
                Vec_DFF_all = DFF.DFF(:,:,final_idx);
                
                DFF_Z = squeeze( ( Vec_DFF_all - nanmean(nanmean(Vec_DFF_all )) )...
                    ./ nanstd(nanstd(Vec_DFF_all)));
                % baseline correction
                
                B_DFF_Z = repmat(nanmean(DFF_Z(1:30,:,:)),[size(DFF_Z,1),1,1]);
                Fluoro  =  DFF_Z - B_DFF_Z ;
                
              clear DFF_Z
              clear Vec_DFF_All
              clear B_DFF_Z
              
            end  
        else
            Fluoro = DFF.DFF(:,:,final_idx);
        end
        
        
        %% plotting 
         neurons = size(Fluoro,3);
         
         plotDFFs(Fluoro,'',neurons)
    end        
end 

end

function plotDFFs(DFF,ClassName,Neurons)

    %% plotting 
figure
hold on 

for nn = 1:size(DFF,3)
    shadedErrorBar([],squeeze(nanmean(DFF(:,:,nn),2))...
        ,squeeze(nanstd(DFF(:,:,nn),[],2))./ sqrt(size(DFF,2))  * 1.96)
end

title(sprintf('%s Responsive: %d Neurons',ClassName,Neurons),...
      'Interpreter','none');
    


end
