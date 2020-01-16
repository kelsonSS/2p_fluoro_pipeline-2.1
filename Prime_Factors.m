function F = Prime_Factors(num)

F = [];

 while mod(num,2) == 0;
     F(end+1) = 2;
     num = num/2;
 end 
 
 
 for ii = 3:2:ceil(sqrt(num))
     while mod(num,ii) == 0
         F(end+1) = ii ;
         num = num/ii  ;
     end 
 end 
         
 if num>2
     F(end+1) = num;
 end 
 
end 
     