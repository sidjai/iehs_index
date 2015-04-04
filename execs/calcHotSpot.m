function [ u ] = calcHotSpot( u, p )
%calcHotSpot Barebones aggregation to get the HotSpot method
%   Inputs is a unit that has gone through the Aggregation
%   outputs the justification output

paraLen = length(u(1).IntVal);

unitLen = size(u,1);

% condEff{1}=cell(condLen,6,paraLen);


for n=1:unitLen
    temp = [];
    condLen = u(n).numCon;
    for cnd = 1:condLen
        nc = sub2ind(size(u),n,cnd);
        [out, tsum] = aggregate('quick',u,nc,p);
        temp= [temp,out{3}];
    end
    HS = abs(mean(temp(:,2:end),2)-temp(:,1));
    HS(logical(isnan(HS)))=0;
    for p=1:paraLen
        k = sort([-.5,HS(p)/2,.5]);
        HS(p)=k(2);
    end
        
    u(n).HotSpot=HS;
end

end

