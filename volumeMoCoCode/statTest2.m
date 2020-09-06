%% General statistical test for comparing two distributions based on normality and distribution size
% May not apply to every scenario -- **Always check the test assumptions**
function [p,resultString] = statTest2(input1,input2)

    if length(input1) == length(input2) % Check if the distributions have same N
        if lillietest(input1) || lillietest(input2) % Check normality
            [p,~] = kruskalwallis([input1 input2],[],'off');
            resultString = 'Non-normal distribution, used Kruskal-Wallis.';
        else
            [p,~] = anova1([input1 input2],[],'off');
            resultString = 'Both distributions were normal, used ANOVA.';
        end
    else % If distributions are different sizes
        if lillietest(input1) || lillietest(input2) % Check normality
            [p,~] = ranksum(input1,input2);
            resultString = 'Non-normal distribution, used Ranksum.';
        else
            [~,p] = ttest2(input1,input2);
            resultString = 'Both distributions were normal, used ttest2.';
        end
    end

end

