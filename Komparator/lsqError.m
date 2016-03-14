function [lsqE] = lsqError(expData,simData)
%function that computes the least squared error from experimental and
%simulated data. sizes of datasets need to fit.
%   Detailed explanation goes here
[nt,nd]=size(expData);
lsqE=0;
expData(isnan(expData))=0;
simData(isnan(simData))=0;
for d=1:nd
    lsqE = lsqE+sum((simData(:,d)-expData(:,d)).^2);
end
end

