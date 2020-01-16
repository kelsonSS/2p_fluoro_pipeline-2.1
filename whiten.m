function [X] = whiten(X)
X=double(X);
X = bsxfun(@minus, X, mean(X));
if  isa(class(X), 'gpuArray');
    
    
    X = X*V*diag(power(complex(1./(diag(D)+eps)),complex((1/2))*V'));
else 


A = X'*X;
[V,D] = eig(A);

    X = X*V*diag(1./(diag(D)+eps).^(1/2))*V';
end