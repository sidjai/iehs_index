function [ s ] = calcUUvals ( u, chem, s, sind, p , j )
%calcUUvals Iterates the index under different upset conditions and
%calculates the UU method
%   UU's are a parameter of the streams due to them being upsetted trying
%   to understand how dangerous that connection is which is a function of
%   the upstream and downstream units. This can be related back to the
%   units using UUstress which is the amount of UU "given" or "absorbed" by
%   the unit.

devVars ={'T','P','F'};
gword = [1 ;-1];
gmean = {'More','Less'};
untLen=size(u,1);
parLen=size(u(1).IntVal,1);
 k=1;
    n=s(sind).to;
    for vind = 1:length(devVars)
        for gind = 1:length(gword)
            if k ~= 1
                uup = struct();
            end
            [uBlock, percent] = makeUpset(u,s,sind,devVars{vind},gword(gind));
            uup = u(:,1);
            if ((length(fields(uBlock)))-(length(fields(u(1)))) ~=0)
                disp('Upset added fields, abandon ship')
            end
            uup(n) = uBlock;
%             for uind = 1:untLen
%                 if uind ==n
%                     uup(n,1)=uBlock;
%                 else
%                     uup(uind,1)=u(uind); 
%                 end
%                 
%             end
            
            
            
            [uup, junk] = calcIntVal(uup,chem,n,1,p,j);
            
            [out, junk] = aggregate('quick',uup,n,p);
            
%             worind=sumup{3,1}; %the mean value ind
            %log of the percent change between the cdps of the upset and
            %normal cases. 
            %percent change between the upset and the original
            change = abs(out{3}(:,1)-u(n).cdp(:,1))./abs(u(n).cdp(:,1));
%             change=abs(cdpup{worind})./abs(u(n).cdp{worind});
            
            
            
            %-ones(parLen,1);
            temp= (change);
            for t=1:length(temp)
                if ~isnan(temp(t))
                    we=sort([-1 temp(t) 1]);
                else
                    we=[0 0 0];
                end
                temp(t)=we(2);
            end
            
            s(sind).UUval{k,1}=temp;
            para= nonzeros(s(sind).UUval{k});
            if isempty(para)
                s(sind).UUavgp(k,1)=0;
                s(sind).UUstdp(k,1)=0;
            else
                s(sind).UUavgp(k,1)=mean(para);
                s(sind).UUstdp(k,1)=std(para);
            end
            
            s(sind).unitUUval{k,1}= gmean{gind};
            s(sind).unitUUval{k,2}= devVars{vind};
            s(sind).unitUUval{k,3}= percent;
            k=k+1;
        end
    end
    UVlen=k-1;
    %Get the average of all the parameters for each of the deviations
    [B,I]=sort(s(sind).UUavgp);
    s(sind).UUwor=[B(end),s(sind).UUstdp(I(end)),I(end)];
    
    %Get the average of all the deviations for each of the parameters
    for p=1:parLen
        temp=0;
        num=UVlen;
        for uv=1:UVlen
            add=chVal(s(sind).UUval{uv}(p),1);
            if add ==0 || isnan(add)
                num=num-1;
            else
                temp= temp + add;
            end
        end
        if num==0
            s(sind).UUs(p,1) = 0; 
        else
            s(sind).UUs(p,1) = temp./num;
        end
    end
end

function val =chVal (varargin)
if length(varargin)>1
    flag=1;
else
    flag=0;
end

in=real(varargin{1});
if flag
    temp = sort([-1,in,1]);
    val= temp(2);
else
    temp = sort([10^-10,in,1]);
    val= temp(2);
    
end
end