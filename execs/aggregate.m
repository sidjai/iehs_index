function [out tsum] = aggregate(varargin)
%aggregate Takes the indivdual methods and makes sense of them
%   First input is a string telling if you want to just take a quick avg or
%   take all the different weights "quick" signifies the quick procedure
%   second is the unit structure
%   third is the streams if you want to aggregate the units together or a
%   unit number, n , to signify the unit you want to aggregate.
%
%   Looks at the userweights provided in unitAlteration worksheet that got
%   parsed
%   Outputs a group of cells that define the different steps in the unit
%   index steps as shown below:
%   1)SecEff add, EDP; 2)Agro chems, EDPAchems; 3)Mix add, MDP
%   4)Vessel add, CDP; 6)Agro units, CDPAunits; 7)UUs add, IDP


u=varargin{2};
W=u(1).userWeights;
chemNames = u(1).chemNames;
paraNames = varargin{4};
if isstruct(varargin{3})
    s=varargin{3};
    flagu=1;

else
    n = varargin{3};
    flagu=0;
end

hierarchy = {'edp'; 'mdp';'cdp';'idp'};
add = {'mix';'Vessel'; 'UUs'};
agroType = {'Average','Best Case','Worst case'};

head={'Mean','std','Median','Lowest','Highest'};

paraLen = length(u(1).IntVal);
unitLen = size(u,1);


%1)SecEff add, EDP; 2)Agro chems, EDPAchems; 3)Mix add, MDP
%4)Vessel add, CDP; 6)Agro units, CDPAunits; 7)UUs add, IDP



%set up the overall "map" of what to do

%condition / chem agro

%do the quick and dirty version first
if ~isempty(strfind(varargin{1},'quick'))

    tmdp = mean(u(n).edp,2)+ u(n).mix;
    tcdp = tmdp + u(n).Vessel;
    out = {u(n).edp, tmdp, tcdp};
    if ~isempty(strfind(varargin{1},'unit'))
        out{4} = tcdp+u(n).HotSpot;
    end
    tsum = 'quick';
else

    if ~flagu
        %single unit aggregation (minus Hotspot)
        for agr=1:3
            tedp = u(n).edp;
            %chemical weights
            [tmdp(:,agr), tsum{agr}] = doWeight(u(n).IntVal,agr,2,paraNames,chemNames);
            tmdp(:,agr) = tmdp(:,agr)+u(n).mix;

            tcdp(:,agr) = tmdp(:,agr)+u(n).Vessel;
        end

        if iscell(W)
            %do the user generated weights
            %1)SecEff add, EDP; 2)Agro chems, EDPAchems; 3)Mix add, MDP
            %4)Vessel add, CDP; 5) add Hotspot, UDP 6)Agro units, UDPAunits; 7)UUs add, IDP

            %redo to parse which weights are used
            tedp = u(n).IntVal+(W{1}.*u(n).fate);
            tmdp(:,4) = W{2}.*(tedp(4)) + W{3}.*u(n).mix;
            tcdp(:,4) = tmdp(:,4)+ W{4}.*u(n).Vessel;
%             out = {tedp, tmdp, tcdp};

        end
        out ={tedp,tmdp,tcdp};

        %add previously calc'ed Hotspot val
        if ~isempty(strfind(varargin{1},'unit'))
            for h=1:size(tcdp,2)
                if h==4
                    out{4}(:,h) = tcdp(:,h)+W{5}.*u(n).HotSpot;
                else
                    out{4}(:,h) = tcdp(:,h)+u(n).HotSpot;

                end

            end
        end



    else
        %unit aggregation
        UUs = [s.UUavgdev];

        for pi=1:paraLen
            vec = UUs(pi,:);
            mUUs(pi,1) = mean(nonZo(vec));
        end


        for n=1:unitLen
            tudp(:,:,n) = u(n).udp;
            unitNames{n}=u(n).name;
        end
        k=1;
        for prev=1:3
            if unitLen > 1
                proMat = tudp(:,prev,:)
            else
                proMat = tudp(:,prev);
            end

            for uagr=1:3


                [tidp{k}, tsum{k}] = doWeight(proMat,uagr,3,paraNames,unitNames);
                tidp{k} = tidp{k} + mUUs;
                tsum{k}{3,1}=['Chem weight: ' agroType{prev} '| Unit weight: ' agroType{uagr}];

                k=k+1;

            end
        end
        if iscell(W)
            for n=1:unitLen
                tudp(:,prev,n) = u(n).cdp(:,4)+W{5}(:,n).*u(n).HotSpot;
            end
            tidp{k} = sum(W{6}.*tcdp(:,4,n),3) + W{7}.*mUUs + u(n).HotSpot;

        end


        out = tidp;

    end

end


end

function [val,smry] = doWeight(mat,agro,dm,stNames,chNames)
staLen=size(mat,1);
if dm > length(size(mat))
  %just return
  val = mat;
  smry{1,1} = ones(staLen,1);
else
  comLen = size(mat,dm);

  switch agro
      case 1
          %average
          val = mean(mat,dm);
      case 2
          %Best case
          [val, smry{1,1}] = min(mat,[],dm);
      case 3
          %worse case
          [val smry{1,1}]= max(mat,[],dm);
  end
end

nval = nonZo(val);

smry{2,1} = [mean(nval), std(nval)];
[smry{2,1}(4), ind] = max(val);
smry{2,2}(4) = stNames(ind);


[smry{2,1}(3), ind] = min(nval);
smry{2,2}(3) = stNames(ind);

% get the names of the values to make sense of them
if agro~=1
    for s=1:staLen

        smry{1,2}{s,1} = chNames{smry{1,1}(s)};
    end %do outside maybe
end

end

function [ weight ] = findWeights( values, type )
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
if type ==1;

    weight = avg/secLen;
else
    weight = worse;
end
end
