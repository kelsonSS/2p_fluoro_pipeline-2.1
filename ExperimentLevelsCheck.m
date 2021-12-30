function Out = ExperimentLevelsCheck(TNBehavior)


 

Out = zeros(length(TNBehavior),2);
    for expt_idx =1:length(TNBehavior)
        
        
        expt= TNBehavior{expt_idx};
        handles = expt.handles{1};
        uLevels = handles.uLevels - 50;
       
        
        
        F = expt.DFF_Z(:,:,expt.Clean_idx & expt.active{:,2}>0);
        
        hits_idx = handles.Hits';
        miss_idx = handles.Miss';
        early_idx = handles.Early';
        
        size(F,2) 
        Out(expt_idx,1) = length(handles.Levels);
        Out(expt_idx,2) = size(F,2);
       
    end 
     Out(:,1) - Out(:,2)