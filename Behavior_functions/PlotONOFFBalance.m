function PlotONOFFBalance(Temporal,SaveName)

if ~exist('SaveName','var')
    SaveName = [];
end 
YoungNumerator =Temporal.YoungGood.ToneBalanceNumber + Temporal.YoungBad.ToneBalanceNumber;

OldNumerator = Temporal.OldGood.ToneBalanceNumber + Temporal.OldBad.ToneBalanceNumber;
  
YoungDenominator = length(Temporal.YoungGood.ToneBiasByAnimal) +...
                   length(Temporal.YoungBad.ToneBiasByAnimal);

OldDenominator = length(Temporal.OldGood.ToneBiasByAnimal) +...
                 length(Temporal.OldBad.ToneBiasByAnimal);
             
YoungPrc = YoungNumerator / YoungDenominator * 100;
OldPrc = OldNumerator / OldDenominator * 100;

Young_label = sprintf('Young: %d / %d', YoungNumerator,YoungDenominator);
Old_label = sprintf('Old: %d / %d', OldNumerator,OldDenominator);

figure
bar(diag([YoungPrc,OldPrc]),'stacked')

title( ' |Tone Bias| < 0.5')
legend(Young_label, Old_label)
xticklabels({'Young','Old'})
ylabel('%')
ylim([0 100])

if SaveName
    saveas(gcf,[SaveName '-TemporalBalance.pdf'])
end 



             
             