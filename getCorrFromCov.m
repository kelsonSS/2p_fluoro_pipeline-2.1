function Output = getCorrFromCov(CovMatrix)
 % This function takes the output from cov()and returns the Psignal or
 % Noise Correlation for each IJ pair of the matrix. 
 
 % 
 % Peusdo-code explination:
 % for all IJ pairs:
 % Corr_ij = cov(ij) / sqrt( var(i) * var(j)) 
 % since covMatrix(i,i) == var(i) && covMatrix(j,j) == var(j)
 % Corr_ij = covMatrix(i,j) / sqrt(covMatrix(i,i) * covMatrix(j,j)) 
 
 
 L = length(CovMatrix);
 
 %pairs = nchoosek(1:L,2)
 
Output = nan(L,L);
for xtemp = 1:L
    for ytemp = xtemp+1:L % this will create a triangular matrix an ignore duplicates
    try
  % xtemp = pairs(ii,1);
  % ytemp = pairs(ii,2);
   % signal correlation equation
   Output(xtemp,ytemp) = CovMatrix(xtemp,ytemp) /...
       sqrt( CovMatrix(xtemp,xtemp) *  CovMatrix(ytemp,ytemp) );
    catch % experiments will be skipped and answer will be nan 
    end 
   end
end 
 
 
 