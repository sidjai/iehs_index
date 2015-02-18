function [ stmTable ] = makeStreamTable( s, p,unitNames )
%makeStreamTable Summuarize stream information into a cell 
%   orders streams and gives the ouput labels
%   adds all the parameters to each stream to show what is going wrong
flag=~estVar(s,'UUwor');

for w=1:length(s)
    if flag(w)
        savgp(w,:)=s(w).UUwor;
    end
end

[B I] =sortrows(savgp,-1);
for sr=1:length(s)
    if sr==1
        stmTable={'Rank','Stream Name','Stream number','Unit from','Unit to',...
            'UU index, worse','Std of other deviations',...
            'Guide word','Variable','Deviation',''};
        for pi = 12: (length(p)+11)
            stmTable{1,pi} = p{pi-11};
        end
        rr=2;
    end
    if sr>length(I)
        stream=sr;
    else
        stream=I(sr);
    end
    
    t=s(stream);
    
    stmTable{rr,1}=sr;
    stmTable{rr,2}=t.name;
    stmTable{rr,3}=stream;
    if t.from==0
        stmTable{rr,4}='Feed';
    else
        stmTable{rr,4}=unitNames{t.from};
    end
    
    if t.to == 999
        stmTable{rr,5}='Product';
    else
        stmTable{rr,5}=unitNames{t.to};
    end
    
    if size(B,1)>=sr && B(sr,1)~=0
        
        stmTable{rr,6}=rund(B(sr,1));
        stmTable{rr,7}=rund(B(sr,2));
        stmTable{rr,8}=t.unitUUval{B(sr,3),1};
        stmTable{rr,9}=t.unitUUval{B(sr,3),2};
        stmTable{rr,10}=rund(t.unitUUval{B(sr,3),3});
        for pi = 12: (length(p)+11)
            stmTable{rr,pi} = rund(t.UUs(pi-11));
        end
    elseif t.to==999
        stmTable{rr,6}='Excluded cause it is a product flow';
    else
        stmTable{rr,6}='Excluded cause Psuedo combined stream';
    end
    rr=rr+1;
        
end

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