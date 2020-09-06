% Compute pairwise absolute distances between points.
%       Inputs: 'X' -- x-coordinates of points
%               'Y' -- y-coordinates of points
%
%       Outputs: 'D' -- matrix of pairwise distances
% ZCB

function D = pairwiseDist(X,Y)

D = NaN(length(X),length(Y));
for i = 1:length(X)
    for j = i:length(Y)
        D(i,j) = hypot((X(i)-X(j)),(Y(i)-Y(j)));
    end
end

end