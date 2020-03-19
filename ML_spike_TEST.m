par = tps_mlspikes('par') % default params
fs = 30;% framerate in hertz
par.dt = .03; % framerate in seconds
par.a = 0.07; % DF/F for one spike
par.tau = 1; % decay time constant (second)



% load fluoro here

% 


% baseline estimation be using moving 30 second average 
wind_size = 30; % window length
wind = fs * wind_size;

baseline_raw  = movmean(Output.fluoAllRaw,fs * wind);
baseline_corr = movmean(Output.fluoAllCorr,fs* wind);

% movemean is not correct for 1:wind so correct with first stable point
baseline_raw(1:wind,:) = repmat(baseline_raw(wind+1,:),[wind,1]);
baseline_corr(1:wind,:) = repmat(baseline_corr(wind+1,:),[wind,1]);


calcium_raw = Output.fluoAllRaw ./ baseline_raw;
calcium_corr = Output.fluoAllCorr ./ baseline_corr;

for ii = 1:10
[spikest{ii} fit drift] = spk_est(calcium_raw(:,ii),par);

[spikest_corr{ii} fit drift] = spk_est(calcium_corr(:,ii),par);
end