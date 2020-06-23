function TN =  GCanalClasses(TN)

for ii = 1:length(TN)
    ClassName  = matlab.lang.makeValidName(TN.Classes{ii});
    classidx = TN.Class_idx == ii;
    classidx = [ii;classidx];
    
    TN.GC.(ClassName) = GCanalTrialsModBalanced_TN(TN,classidx);
end 