function [AnimalInfo] = AnimalAgeAnalysis(AnimalInfo,StartDate,EndDate)
%%% this function takes in a AnimalInfo object and the starting and
%%% ending Dates extracted from the Psigal files from a set of experiments 
%%% and returns the age of the animal at the start and end of that set of
%%% experiments. if only one Date is given the same file will be 
%%% used for the start and end dates


if ~exist('EndDate','var')
    warning(sprintf(['no EndDate given:',... 
        '\n using StartDate as  Start- and  EndDate']))
    EndDate = StartDate;
end 



     AnimalInfo.StartAgeMonths = caldiff([AnimalInfo.DOB,...
                                   StartDate] ,'Months') ;
     AnimalInfo.EndAgeMonths = caldiff([AnimalInfo.DOB,...
                                   EndDate] ,'Months') ;
    
     AnimalInfo.StartAgeExact = StartDate - AnimalInfo.DOB;
     AnimalInfo.EndAgeExact =  EndDate - AnimalInfo.DOB;
     AnimalInfo.MeanAgeDays = days(mean([AnimalInfo.StartAgeExact,...
                                         AnimalInfo.EndAgeExact]));
    