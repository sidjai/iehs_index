function [ tab ] = makeUnitTable( u, p)
%makeUnitTable Orders the Units by decreasing UDP from left to right
%   Also includes the average of each of the methods below the value along
%   with some statistical analysis

for n=1:length(u)
    udps(n,1)= mean(u(n).udp(:,1));
%         savgp(w,:)=s(w).UUwor;
    
end

[B I] =sort(udps,'descend');
% untTable = {'Unit Number','Unit name','Unit Type',...
%     'CDP, Avg','CDP, std', 'MDP, avg', 'MDP, std','EDP','IntVal','Highest Hazard'};

labels ={'Number','Name','Type','UDPavg parameters and aggregation','UDPmaxParameter','UDPstd of parameters'...
    '','Chem','SecEffs','Mix','Vessel','HotSpot',...
    '','ChemPara','Chemstd','SecEffsPara','SecEffsStd','MixPara','MixStd','VesselPara','VesselStd','HotSpotPara','HotSpotStd'};
for c=1:length(u)+1
    
    if c==1
        for ab=1:length(labels)
            tab{ab,1} = labels{ab};
        end
    else
        n=I(c-1);
        t=u(n,1);
        tab{1,c} = n;
        tab{2,c} = t.name;
        tab{3,c} = t.type;
        [tab{4,c}, tab{5,c} tab{6,c}] = getStats(t.udp(:,1),p);
        
        [tab{8,c}, tab{14,c} tab{15,c}] = getStats(t.IntValRaw,p);
        [tab{9,c}, tab{16,c} tab{17,c}] = getStats(t.fateRaw,p);
        [tab{10,c}, tab{18,c} tab{19,c}] = getStats(t.mix,p);
        [tab{11,c}, tab{20,c} tab{21,c}] = getStats(t.Vessel,p);
        [tab{12,c}, tab{22,c} tab{23,c}] = getStats(t.HotSpot,p);
        
    end
    
           
    
%     untTable{ind,1}=n;
%     untTable{ind,2}=t.name;
%     untTable{ind,3}=t.type;
%     untTable{ind,4}=rund(B(ind-1,1));
%     untTable{ind,5}=rund(B(ind-1,2));
%     untTable{ind,6}=rund(t.summdp{2,1});
%     untTable{ind,7}=rund(t.summdp{2,2});
%     edp=nonzeros(t.edpRaw);
%     untTable{ind,8}=rund(mean(edp));
%     ini=nonzeros(t.IntValRaw);
%     untTable{ind,9}=rund(mean(ini));
%     
%     [junk para]=max(t.cdp{B(ind-1,3)});
%     untTable{ind,10}=p{para};
    
        
end

end
function [nu, para, sig] = getStats(vec,p)
[junk, k] = max(vec);
para = p{k};
nu = rund(mean(vec));
sig = rund(std(vec));
end

function out = rund(in)
messy=floor(in*1000)/1000;
if messy==0
    if in==0
        out=0;
    else
        out=10^-10;
    end
else
    str=num2str(in);
    ow=textscan(str,'%5.3f');
    out=ow{1}(1);
end
% if out==0 && in~=0
%     out=10^-10;
% end
end