% Compute pairwise absolute distances between points.
%       Inputs: 'X' -- x-coordinates of points
%               'Y' -- y-coordinates of points
%
%       Outputs: 'A' -- matrix of pairwise angles (degrees)
% ZCB

function A = pairwiseAngle(X,Y)

A = NaN(length(X),length(Y));
for i = 1:length(X)
    for j = i+1:length(Y)
        A(i,j) = atan2d((X(i)-X(j)),(Y(i)-Y(j)));
    end
end

end