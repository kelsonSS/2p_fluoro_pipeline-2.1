function DF =  Munge_DF(DF,norm_mode)

    if ~exist('norm_mode','var')
        norm_mode = 'Z-score';
    end 


if iscell(DF) % if Cell iterate over cells to get full index 
    temp = DF;
    clear  DF
    DF = struct()
    DF.DFF_norm = []; 
  %  DF.DFF_norm_mu= [];
    DF.Clean_idx = [];
    DF.active = [];
    DF.experiment_list = [];
    for ii =1:length(temp)
        expt = Munge_DF(temp{ii},norm_mode );
         % append fields to master list 
        
         if ii> 1 && size(expt,2) ~= size(DF.DFF_norm,2);
             if ~iscell(DF.DFF_norm)
                 DF.DFF_norm = {DF.DFF_norm};
               %  DF.DFF_norm_mu = {DF.DFF_norm_mu};
             end 
             DF.DFF_norm{ii} = expt.DFF_norm;
%              DF.DFF_norm_mu{ii} = expt.DFF_norm_mu;
         else
            DF.DFF_norm = cat(2,DF.DFF_norm, expt.DFF_norm  );
          % DF.DFF_norm_mu = cat(2,DF.DFF_norm_mu, expt.DFF_norm_mu  );
            
         end 
        DF.Clean_idx = cat(1,DF.Clean_idx,temp{ii}.Clean_idx);
        DF.active = cat(1,DF.active,temp{ii}.active);
        DF.experiment_list= cat(1,DF.experiment_list,...
                           ones(size(expt.DFF_norm,3),1) * ii);
    end 
 clear temp 

else 
if isstruct(DF)
    DF.DFF2 = squeeze(nanmean(DF.DFF,2));
else
    temp = DF;
    clear DF
    DF = struct('DFF',temp, 'DFF2', squeeze(nanmean(temp,2)) )
    clear temp
end

 %%  normalizing to absolute max 

 switch norm_mode
    case 'normalized'
        DFF_ab_max = squeeze(max(abs(DF.DFF2)));
         %DF.DFF_norm = DF.DFF_norm./DFF_ab_max;
         DF.DFF_norm = DF.DFF2./DFF_ab_max;

%% normalizing to baseline corrected z-score (decide which one)
    case 'Z-score'
        %DF.DFF_norm = DF.DFF_Z;
        DF.DFF_norm = squeeze(nanmean(DF.DFF_Z,2));
      
    otherwise
        error('norm_mode must be set to normalized or Z-score')
 end 
end 
