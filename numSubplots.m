function [ rc ] = numSubplots(n)
   % this function takes in the number of total plots to display and 
   % gives back a two element matrix consisting of the correct number 
   % of rows and columns to use in the subplot function 
   % 
   
   % ex:  rc = numSubplots(6)
       %   rc = [2 3]
       %    figure
       %     subplot(rc(1), rc(2),1)
       %     ... %plot one here 
       %
       
  if ~ isnumeric(n) ||( floor(n) ~= (n) || n < 0)     
       errordlg('Must imput an integer number')
  end 
  
  if n < 4
      rc = [1 n];
      return
  end 
  
   x = sqrt(n);
    
   if  x-floor(x) < .5 
    rc(1) = floor(x); 
    rc(2) = ceil(x);
   else 
    rc(1:2) = ceil(x);   
   end 
   
end 