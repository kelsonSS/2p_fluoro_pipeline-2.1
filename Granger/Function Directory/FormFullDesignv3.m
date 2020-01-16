function [Xneff,Dnf] = FormFullDesignv3(Xcross,Xself,ct) 
%%% v3: Trial-based Full design matrix

[Np,Mhc,Ncells,R] = size(Xcross);
cellC = 1:Ncells; cellC(ct) = [];

%%% FULL Model
% Form Full model design/covariate matrix Xf
Xfc = Xself;
for cc = 1:Ncells-1
    Xfc = cat(2,Xfc,squeeze(Xcross(:,:,cellC(cc),:)));
end
Xf = Xfc; 
%%% FORM an EFFECTIVE covariate matrix
% Standardize the effective covariate matrix Xeff
Mf = size(Xf,2);
Xeff = zeros(R*Np,Mf); 
for r = 1:R
    Xeff((r-1)*Np+1:r*Np,:) = Xf(:,:,r);
end  % can use cat(1,...) instead
Dnf = diag(sqrt(var(Xeff)));
Xneff = Xeff/Dnf;
   