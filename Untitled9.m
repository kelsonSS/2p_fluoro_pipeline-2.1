[ rc ] = function numSubplots(n)
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
       
  if ~isinteger(n)     
       errordlg('Must imput an integer number')
  end 
  
   x = sqrt(n);
    
   rc(1) = floor(x); 
   rc(2) = ciel(x);
       