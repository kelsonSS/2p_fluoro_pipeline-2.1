function C = CreateTreeClassifier(Neurons,FreqLevelOrder,type )


if ~exist('type','var')
    type = 'Freqs'
end 
 
switch type 
    case 'Freqs'
         FreqLevelOrder = FreqLevelOrder(:,1); 
    case 'Levels'
         FreqLevelOrder = FreqLevelOrder(:,2);
    case 'TonesInQuiet'
         type = 'Freqs';
         
         Idx = FreqLevelOrder{:,2} == 99
         FreqLevelOrder = FreqLevelOrder(Idx,1);
         Neurons = Neurons(Idx,:);
        
    otherwise 
        return 
end 

 F =  [array2table(Neurons),  FreqLevelOrder];
 
 t = templateTree('surrogate','on','MinLeafsize',9);
 
 C = fitcensemble(F,type,'OptimizeHyperParameters','auto')
 
 %C = fitcensemble(F,type,'Method','AdaBoostM2',...
 %    'NumLearningCycles',400,...
 %    'LearnRate',.18196,...
 %    'Learners',t,...
 %    'CrossVal','on')  ;

kflc = kfoldLoss(C,'Mode','cumulative');
figure;
plot(kflc);
ylabel('10-fold Misclassification rate');
xlabel('Learning cycle');
 
end 