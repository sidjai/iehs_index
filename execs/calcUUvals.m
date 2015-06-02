function [ s ] = calcUUvals ( u, chem, s, sind, p , j )
%calcUUvals Iterates the index under different upset conditions and
%calculates the UU method
%   UU's are a parameter of the streams due to them being upsetted trying
%   to understand how dangerous that connection is which is a function of
%   the upstream and downstream units. This can be related back to the
%   units using UUstress which is the amount of UU "given" or "absorbed" by
%   the unit.

devVars ={'T','P','F'};
gword = [1; -1];
gmean = {'More','Less'};
untLen=size(u,1);
parLen=length(p);
k=1;
upsetTotTime = 15;

for vind = 1:length(devVars)
    for gind = 1:length(gword)
        if k ~= 1
            uup = struct();
        end
        upsetTime = upsetTotTime; %min
        uup = u(:,1);

        upsetSet = [];
        upsetSetRes = [];
        orgRes = [];
        perSet = [];

        while upsetTime>0
%             if upsetTime<15
%                 upsetTime = upsetTime;
%             end
            if upsetTime==upsetTotTime
                sdev = s(sind);
                upType = gword(gind);
            else
                sdev = s(u(badU).out(1));
                upType = badUnit2badStream(u(badU), sdev);

%                 sdev = s(u(badU).out(1));

            end

            badU = sdev.to;

            %break if the upsets made a loop to the same unit
            %break if you hit a product flow cause there isn't a unit to
            %upset
            if (ismember(badU,upsetSet)) || badU==999
                break;
            end

            upsetSet = [upsetSet badU];

            indvUpsetTime = min(u(badU).tau,upsetTime);

            [uBlock, percent] = makeUpset(uup,sdev,...
                devVars{vind},upType,indvUpsetTime);
            perSet = [perSet percent];
            
            if ((length(fieldnames(uBlock)))-(length(fieldnames(u(1)))) ~=0)
                disp('Upset added fields, abandon ship')
            end

            uup(badU) = uBlock;

            [uup, junk] = calcIntVal(uup,chem,badU,1,p,j);
            [out, junk] = aggregate('quick',uup,badU,p);
            upsetSetRes = [upsetSetRes out{3}(:,1)];
            orgRes = [orgRes u(badU).cdp(:,1)];

            upsetTime = upsetTime - indvUpsetTime;

        end



        %percent change between the upset and the original
        change = (upsetSetRes-orgRes);
%         change = abs(out{3}(:,1)-u(n).cdp(:,1))./abs(u(n).cdp(:,1));


        change = chVal(change);

        %Add up all the different units that were effected
        %Use sum cause this represents the total effect of 1 upset
        %(para x units) -> (para x 1)
        allUCh = sum(change,2);

        s(sind).UUval{1,k} = change;
        s(sind).UUsummult(:,k) = allUCh;

        s(sind).unitUUval{1,k} = gmean{gind};
        s(sind).unitUUval{2,k} = devVars{vind};
        s(sind).unitUUval{3,k} = perSet;
        s(sind).unitUUval{4,k} = upsetSet;
        k=k+1;
    end
end

UUmat = s(sind).UUsummult;

%Get the average of all the parameters for each of the deviations
%dim(1 x deviations)
for di = 1:size(UUmat,2)
    s(sind).UUavgp(1,di) = mean(nonZo(UUmat(:,di)));
    s(sind).UUstdp(1,di) = std(nonZo(UUmat(:,di)),0,1);
    [B,I]=sort(s(sind).UUavgp);
    s(sind).UUwor=[B(end),s(sind).UUstdp(I(end)),I(end)];
end

%Get the average of all the deviations for each of the parameters
%dim(para x1)
for pi = 1:parLen
    s(sind).UUavgdev(pi,1) = mean(nonZo(UUmat(pi,:)));
    s(sind).UUstddev(pi,1) = std(nonZo(UUmat(pi,:)),0,2);
end

end
function sout = badUnit2badStream(uin, sorg)

sout = sorg;

chFields = {'CP','T','rho','x','H','P'};

for hf = 1:length(chFields)
    sout.(chFields{hf}) = uin.(chFields{hf});
end
end

function val =chVal (varargin)
if length(varargin)>1
    lowThres = -1;
    lowVal = lowThres;
else
    lowThres = 0;
    lowVal = 10^-10;
end

val=real(varargin{1});
nutsSet = ~isfinite(val);
overSet = (val>1);
underSet = (val < lowThres);
val(nutsSet) = 0;
val(overSet) = 1;
val(underSet) = lowVal;

end
