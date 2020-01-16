 corrs = function getCorrFromCov(CovMatrix)
 % This function takes the output from cov()and returns the Psignal or
 % Noise Correlation for each IJ pair of the matrix. 
 % 
 % Peusdo-code explination:
 % for all IJ pairs:
 % Corr_ij = cov(ij) / sqrt( var(i) * var(j)) 
 % since covMatrix(i,i) == var(i) && covMatrix(j,j) == var(j)
 % Corr_ij = covMatrix(i,j) / sqrt(covMatrix(i,i) * covMatrix(j,j)) 
 
 
 
 
 pairs = nchoosek(1:length(CovMatrix),2); 

s_corr = zeros(RANs);
for ii = 1:length(pairs)
   xtemp = pairs(ii,1);
   ytemp = pairs(ii,2);
   % signal correlation equation
   s_corr(xtemp,ytemp) = CovMatrix(xtemp,ytemp) /...
       sqrt( CovMatrix(xtemp,xtemp) *  CovMatrix(ytemp,ytemp) );
 
end 
 
 
 
 
 