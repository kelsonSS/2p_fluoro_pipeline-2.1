function spikes = deconSpikesFromCal(DFF,t,T,sig,lam,tau)
[nroi npoints] = size(DFF);
mtime = max(t);
% simulate data
% T       = 10000; % # of time steps
% N = poissrnd(P.lam*V.dt*ones(T,1)); % simulate spike train
% C = filter(1,[1 -P.gam],N);         % calcium concentration
% F = P.a*C+P.b + P.sig*randn(T,1);   % observations
times = [];
vals = [];
phat=[];
psig = [];
chans = [];
freq = zeros(1, nroi);
if ~isempty(lam)
    P.lam = lam;
end
V.dt = T;
P=[];
rejcell=0;
for j = 1 : nroi
    F = DFF(j,:);
    if ~isempty(tau)
        P.gam = 1-V.dt/tau(j);
    end
    if ~isempty(sig)
        psig = [psig,sig(j)];
        P.sig = sig(j);
    end
    
    %get Nhat from fast_oopsi
    if sum(F)==0
        rejcell=rejcell+1;
        warning('off','MATLAB:nearlySingularMatrix')
        warning('off','MATLAB:illConditionedMatrix')
        disp([num2str(rejcell) ' rejected cell(s): ' num2str(j)])
    end
    close all
    [Nhat, Phat, Vtmp, C] = fast_oopsi(F,V,P);
    phat = [phat, Phat];
    %Find spike indices
    tpos = find(Nhat);
    tpos2 = find(tpos <= numel(t));
    %Spike times
    times = [times,t(tpos(tpos2))];
    %Continuous values for inferred spikes (oopsi does not contrain spike
    %values to be integers)
    vals = [vals, Nhat(tpos(tpos2))'];
    tchans = ones(1,numel(tpos(tpos2)))*j;
    chans = [chans,tchans];
    freq(j) = sum(Nhat(tpos(tpos2)))/mtime;
    if sum(F)==0
        warning('on','MATLAB:nearlySingularMatrix')
        warning('on','MATLAB:illConditionedMatrix')
    end
end
%remove negative data points
tpos = find(vals <= 0);
if ~isempty(tpos)
    warndlg('OOPSI inferred negative spikes!!! oopsidaisy!')
end
vals(tpos) = [];
times(tpos) = [];
chans(tpos) = [];
%Output
spikes.times = times;
spikes.vals = vals;
spikes.phat=phat;
spikes.chans = chans;
spikes.freq = freq;


