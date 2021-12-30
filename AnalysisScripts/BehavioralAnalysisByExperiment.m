function  [resultsAll] =  BehavioralAnalysisByExperiment(ActiveFolders,toPlot)

% this function takes in the AnimalID which corresponds to a folder
% containing the psignalfiles of this animals behavioral performance. This
% function then will look at the daily and overall performance of the
% animal
% Inputs
% AnimalID - name of animal folder containing psignalfiles
% toPlot - optional limit number of figures to plot 
%          All - default
%          Totals - only the combined plots of animals overall 
%                  performance will be plotted 
%
SNRs = [20,10,0]
to_save = 0;
toPlot = 0;
daily_plot_flag = 0;

% TODO : extend to be able to graph data generated from training box
Main_Psignal_Folder = 'Z:\Kelson\Analyzed\';
ExptResults = {}
psigAll = {}
for expt_idx = 1:length(ActiveFolders)  
    % get psignal info
    psig_file = ActiveFolders{expt_idx}
       
     
    psigDay = WF_getPsignalInfo(psig_file );
   
    % check for correct sound class
    % if  ~contains(psigDay.Class,'Tone')
    %     continue
    % end 
    % plot 
    if daily_plot_flag
        figure
        title(PsignalFiles(expt_idx), 'Interpreter','None')
    end 
    if isfield(psigDay,'Performance') % only look at behavior expts.
        expt = ExtractAnimalBehaviorDaily(psigDay,daily_plot_flag);
        % munge into correct form for concatenation
        results.DprimeTotal = expt.DprimeTotal;
        results.Dprime2 = expt.Dprime2;
        results.numTrials = sum(expt.TrialsPerLevel);
        results.HitRate = expt.HitRate(end);
        results.EarlyRate = expt.EarlyRate(end);
        results.EarlyRate2 = expt.EarlyRate2;
        results.MissRate = expt.MissRate(end);
        results.SNRs = {expt.SNR};
        results.PercentCorrectLevel = {expt.PercentCorrect};
        results.dPrimeLevel = {expt.DprimeLevel};

        % Package into Output Structure
        if expt_idx == 1
            resultsAll = results;
        else
            resultsAll = ConcatenateStructs(resultsAll,results,1);
        end
        
    end
end
    

end 
