%TonesNoiseActive Analysis Script



% bandwith analysis 


%Active GC
GCanalTrialsModBalanced_TN(TonesNoiseActive(:,2),'Active')
GCanalTrialsModBalanced_TN(TonesNoiseActive(:,2),'HitMiss')

%Passive GC
GCanalTrialsModBalanced_TN(TonesNoiseActive(:,1),'Passive')


% Correlations

FRA_old_noise.Corr = Correlations(FRA_old_noise);
FRA_young_noise.Corrs = Correlations(FRA_young_noise);

% bayes 
 Bayes = BayesClassifierPassive

