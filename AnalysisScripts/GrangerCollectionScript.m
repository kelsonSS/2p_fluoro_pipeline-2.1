% Granger Results Collection Script 
 
% this script uses GrangerData Object and CollectResults function to
% iteratively collect all the data produced by the Granger Causality code 
% Kelson Shilling-Scrivo

%init
Classes = {'','_Noise','_Tones','_Offset';... % file suffix
           'All','Noise','Tones','Offset'};  % Class names 
FilesTemp = [];
GrangerResults = struct();

% main loop
for class = 1:size(Classes,2)
FilesTemp = cellfun(@(x) [x Classes{1,class}],GrangerData,'UniformOutput',0);

GrangerResults.(Classes{2,class}) = CollectGCResults(FilesTemp)  ;
          
end    
          
          