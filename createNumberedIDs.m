function out = createNumberedIDs(high,reps)
% creates a cell of numbers 1:high repeated 1 x reps times

out =  repmat([1:high],reps,1);

out =  arrayfun(@num2str,out(:),'UniformOutput',0)


