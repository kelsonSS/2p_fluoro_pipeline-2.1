
function [p,resultString] = findSignficance(input1,input2,MultCompare)

% This function extends statTest2 by giving it the ability to do multiple
% comparision tests of the Anova and KruskalWallace test 

if ~exist('MultCompare','var')
    MultCompare = true;
end 

if exist('input2','var')
   p,resultString = statTest2(input1,input2);
end

if length(size(input1))~= 2 
    errordlg('Must input an MxN Matrix')
end 
 
 rc = size(input1); % row-columns
 
 % Analysis functions take a [measurement X condition]  Matrix
 % there should always be more measurements than conditions so if we see
 % the opposite the matrix is likely in the wrong orientation
     if rc(1) < rc(2)  
      input1 = input1';
      rc = size(input1);
     end 
    
     
   normality = 1;     
 for col_idx = 1: rc(2)
     
     if lillietest(input1(:,col_idx))
         normality = 0 ;
         break
     end 
 end 
     
 
 
 
        if ~normality % non normal
            [p,~,stat] = kruskalwallis([input1],[],'off');
            resultString = 'Non-normal distribution, used Kruskal-Wallis.';
        else % normal
            [p,~,stat] = anova1([input1],[],'off');
            resultString = 'Both distributions were normal, used ANOVA.';
        end
        
        if MultCompare
           p = multcompare(stat);
        
        end 
  
end

