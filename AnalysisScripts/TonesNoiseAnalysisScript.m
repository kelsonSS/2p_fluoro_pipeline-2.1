%TonesNoiseActive Analysis Script



% bandwith analysis 


%Active GC
GCanalTrialsModBalanced_TN(TonesNoiseActive(:,2),'Active')
GCanalTrialsModBalanced_TN(TonesNoiseActive(:,2),'HitMiss')

%Passive GC
GCanalTrialsModBalanced_TN(TonesNoiseActive(:,1),'Passive')


% bayes 
 Bayes = BayesClassifierPassive

