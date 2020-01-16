function [wn,sig2] = StaticEstimCV(Xn,robs,gamma,LL,LR)

%%% Static Estimation of Sparse Tuning Coefficients (including Cross-History)
%%% Feb 20th 2017

% schid: Scheme ID

[N,M] = size(Xn);
%muhat = zeros(LR,1);

%%% Estimate the variance using ML estimate
XXn = Xn'*Xn;

Xdn = Xn'*robs;
wML = XXn\Xdn;
sig2ML = (norm(robs - Xn*wML)^2)/N;
% Set Initials based on ML estimate
% muML = mean(robs);
% zx = sum(Xn)';
% wML = XXn\(Xdn-muML*zx);
% sig2ML = (norm(robs -muML - Xn*wML)^2)/N;

% Set Initials based on ML estimate
sig2l = sig2ML;
wl = wML;

for ll = 1:LR
    XXl = XXn/sqrt(sig2l*N*log(M));
    Xdl = Xdn/sqrt(sig2l*N*log(M));     
    % Step-size Selection
    rhom = max(abs(eig(XXl))); al = 0.9/rhom;
    for l = 1:LL
        gl = (Xdl - XXl*wl);
        wl = SoftThreshold(wl + al*gl ,gamma*al); 
    end
    % Compute Variance
    sig2l = (norm(robs - Xn*wl)^2)/N;
    %Sig2ft(ll,ct) = sig2fl;            
end
wn = wl; sig2 = sig2l;

        
end









