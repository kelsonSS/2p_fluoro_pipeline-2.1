function Xz = FormHistMatrixv2(x,Whc,Lm)
% Mod: Mar 27 2017

[N,Ncells] = size(x);
Np = N - Lm;
Mhc = numel(Whc); Lhc = sum(Whc);

X = zeros(Np,Mhc,Ncells);
Lhself1 = zeros(Mhc,1); Lhself2 = Lhself1;
for j = 1:Mhc
    Lhself1(j) = sum(Whc(Mhc-j+2:Mhc));
    Lhself2(j) = sum(Whc(Mhc-j+1:Mhc));
end
% Response History Components WINDOW-BASED
for k = 1:Np
    for cc = 1:Ncells
        for j = 1:Mhc
            X(k,Mhc-j+1,cc) = mean(x(k+Lhself1(j)+Lm-Lhc:k+Lhself2(j)+Lm-Lhc-1,cc));
        end 
    end
end
% Zero-mean History Components
Xz = X;
for cc = 1:Ncells
    Xz(:,:,cc) = X(:,:,cc) - ones(Np,1)*mean(X(:,:,cc)); % Corrected mean on April 4th for window-based history
end


end
