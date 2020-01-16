
function S = recursive_3_sum(range,target,list)

if length(list) == 3 
    if sum(list) == target
        S(end+1,:) = list;
    end
    return 
end 

for ii = 1:length(range)
  trange = range;
    trange(ii) = [];
   S = recursive_3_sum(trange,target,[list range(ii)]);
end 

end 
