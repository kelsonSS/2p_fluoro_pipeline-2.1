function H = myEntropy(s)
% information entropy equation 

% sum probabilities should equal 1
s = s./sum(s(:));

% entropy eq
H = -sum(s .* log2(s)) ; 


