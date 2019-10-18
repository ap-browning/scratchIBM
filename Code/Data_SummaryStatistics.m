function [N,PC,D] = Data_SummaryStatistics(X)
%SummaryStatistics
% Return summary statistics, given agent locations X

    N   = size(X,1);
    PC  = Data_PairCorrelationFcn(X,400,5,100);
    D   = histcounts(X(:,2),0:1900/80:1900);

end