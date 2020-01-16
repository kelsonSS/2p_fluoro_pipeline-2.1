function [wn,sig2,Dev,B] = StaticEstim(Xn,robs,gamma,LL,LR)

%%% Static Estimation of Sparse Tuning Coefficients (including Cross-History)
[N,M] = size(Xn);
%%% Estimate the variance using ML estimate
XXn = Xn'*Xn;
Xdn = Xn'*robs;
% Set Initials based on ML estimate
wML = XXn\Xdn;
sig2ML = (norm(robs - Xn*wML)^2)/N;
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
end
wn = wl; sig2 = sig2l;
% Compute gradient at updated wl and sig2l
gn = (Xdn - XXn*wn)/sig2; 
Hn = XXn/sig2;

% Compute the Unbiased Deviance
B = gn'*(Hn\gn); 
Dev = -N*log(sig2) + B; 
        
end









