function data = clip(data,xmin,xmax)
% this function creates bounds for your data such that for each entry 
% xmin <= data(ii) <= xmax 

if isempty(data)
    return
end 

assert(xmax > xmin)



data(data < xmin) = xmin;
data(data > xmax) = xmax;