function [BFs,BF_Field]= BestFrequencyAnalysisAll(passive)

% extends BestFrequencyAnalysis

for expt = 1:length(passive)

[BFs(expt),BF_Field(expt)] = BestFrequencyAnalysis(passive{expt});

end

