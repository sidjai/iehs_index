function [ weight ] = findWeights( values )
%From a matrix find the average and worst case sceanario weights as a cell
%   Detailed explanation goes here

parLen = size(values,1);
secLen = size(values,2);
avg=ones(parLen,secLen);
worse=zeros(parLen,secLen);
    [temp, worInd] = max(values,[],2);
    for p=1:parLen
        worse(p,worInd(p))=1;
    end

    weight={avg/secLen,worse};
end

