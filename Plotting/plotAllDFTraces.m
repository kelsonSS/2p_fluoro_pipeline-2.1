function plotAllDFTraces(Data)


DFF = Data.DFF;

figure
for ii = 1:1000
subplot(10,10,mod(ii-1,100)+1)
shadedErrorBar([], nanmean(DFF(:,:,ii),2) , nanstd(DFF(:,:,ii),[],2) )
axis off 
if mod(ii,100) == 0
pause
end

end