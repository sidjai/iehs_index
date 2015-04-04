function [ chem ] = loadSafetyInfo( exfile, crw, chem)
%loadSafetyInfo Grabs info from the chemical safety information excel input
%sheet and turns it into the a structure with the given fields
%Also grabs from the analysis report for the chemicals in aspen property
%   1 chem per sheet

% fields = {'Name';'formula';'CAS';'SMILES';'MW';'CPGConst';'CPLConst';'rho';'rhov';'Vmol';'Tb';'Tc';'Tmelt';'Ktherm';'dHvCon';'acentric';'RG';'DOT';'Hcode';'Scode';'antConst';'dHburn';'Flash';'LFL';'UFL';'LOC';'TLVc';'IDLH';'LCInhal';'AIT';'dHdecomp';'dHform';'dGform';'Kst';'overPres';'TESTIrr';'pH';'Ldderm';'TLV';'REL';'TESTChr';'LDChr';'ERPG2';'BCFreal';'BCFpred';'Kow';'Percistency';'Halflifereal';'Halflifepred';'OCED';'BODratio'};

% indVal=4;%future find these indexes
findLine = @(c , string) (find(~cellfun('isempty', strfind(c,string))));

chmLen=length(chem);

%the RG labels.
noncombust = [2;44;100;37;98;38;39];
oxidizers = [2;44];
organics = [33;47;17;27;32;42;20;30];

estFields = {'Tc','Pc','acentric','CPGConst','antConst','dHvCon','dHform','dGform'};
% Load the excel safety sheets and sort them
for c=1:chmLen
    
      [junk,junk,raw]=xlsread(exfile,c);
      
if c==1    
    %find where they put the all the headings
    for h=1:size(raw,2)
        if (ischar(raw{1,h}) || ~isnan(raw{1,h}))
            if ~isempty(strfind(raw{1,h},'ble'))
                indVar=h;
            elseif ~isempty(strfind(raw{1,h},'alue'))
                indVal=h;
            elseif ~isempty(strfind(raw{1,h},'nit'))
                indUnit=h;
            elseif ~isempty(strfind(raw{1,h},'ield'))
                indField=h;
            end
        end
    end
    fnd=1;
    last=size(raw,1);
    for ih=2:last
        if (ischar(raw{ih,indField}) || ~isnan(raw{ih,indField}))
            fields{fnd}=raw{ih,indField};
            fnd=fnd+1;
        end
    end
end
name=raw(3,indVal);

creal=find(strcmp({chem.ID},name{1}));

if ~isempty(creal)
    chem(creal).safetyInfo=raw;
end


end


% Go through the safety excel sheet
for c=1:chmLen
    
  
%     indVal=4;%future find these indexes
    raw=chem(c).safetyInfo;
    
    
    %find the indices with more than one item going horizontally
    ind=2;
    exclude=[];
    while ind<=last
        if ~isnan(raw{ind,indVal+1}) &&  isempty(findLine(estFields,fields{ind-2}))
            exclude=[exclude,ind];
        end
        ind=ind+1;
    end
    %read horizontally for these ones to grab all the data.
    for items=exclude
        cons=1;
        while ~isnan(raw{items,indVal-1+cons})
            chem(c).(fields{items-2})(cons,1)=raw{items,indVal-1+cons};
            cons=cons+1;
        end
    end
    
    %Find the ones in farenheit
    ind=2;
    unitch=[];
    
    while ind<=last
        if (ischar(raw{ind,indUnit}) || ~isnan(raw{ind,indUnit}))
            
            if strcmpi(raw{ind,indUnit},'f')
            raw{ind,indVal}=(raw{ind,indVal}-32)/1.8;
            end
            
            if strcmpi(raw{ind,indUnit},'c') && strcmp(fields{ind-2},'Tc')
                raw{ind,indVal}=(raw{ind,indVal}+273.15);
            end
        end
        
            
        ind=ind+1;
    end
    
    %Equate the safety information provided with the given field.
    ind=2;
    while ind<=last
        temp=raw{ind,indVal};
        if (ischar(temp) || ~isnan(temp)) ...
                && isempty(find(exclude ==ind,1))...
                && (~isfield(chem(c),fields{ind-2}) || isempty(chem(c).(fields{ind-2})))
%                 && isempty(findLine(estFields,fields{ind-2}))

            chem(c).(fields{ind-2})=temp;
        end
        ind=ind+1;
    end
    
    %
    %Additional parameters
    %
    kw=1;
    for r=chem(c).RG'
        teCp(1,kw)=isempty(find(noncombust==r,1));
        teCp(2,kw)=~isempty(find(oxidizers==r,1));            
        teCp(3,kw)=~isempty(find(organics==r,1));
        
        kw=kw+1;
    end
    [ro, co]=find(teCp);
    
    chem(c).flagFuel=~isempty(find(ro==1,1));
    chem(c).flagOx=~isempty(find(ro==2,1));
    chem(c).flagOrg=~isempty(find(ro==3,1));
    
    
    %find the number of hydrogens, etc.
    m=locAtom(chem(c).formula,'C');
    x=locAtom(chem(c).formula,'H');
    y=locAtom(chem(c).formula,'O');
    
    su=locAtom(chem(c).formula,'S');
    chem(c).compStruc(1)=m;
    chem(c).compStruc(2)=x;
    chem(c).compStruc(3)=y;
    chem(c).oxyStoic=m+(x/4)-(y/2);
    chem(c).oxyCont=-1600*(2*m + x/2 -y)/chem(c).MW(1);
    chem(c).contCHNO=isempty(find([m,x,y,locAtom(chem(c).formula,'N')]==0,1));
    
    %rel bouyancy
    if isfield(chem(c),'rho') && ~isempty(chem(c).rho)
        chem(c).bouy=(9.8*(chem(c).rho-1.18)/1.18);
        
    else
        chem(c).rho=chem(c).rhov;
        chem(c).bouy=(9.8*(chem(c).rhov-1));
    end
    %heat of combustion
    
    chem(c).dHburn=(((x/2)*-241818) - (m*393509) - (su*296830)...
        - (chem(c).dHform));%Kj/kmol
    
    %optimum fuel to oxidizer ratio mass
    chem(c).fueloxopt= chem(c).MW/chem(c).oxyStoic*32;
    
    % make CP and Psat into functions
    chem(c).Psat = @(T) (chem(c).antConst(1) + chem(c).antConst(2)/(T+chem(c).antConst(3)) + chem(c).antConst(4)*T + chem(c).antConst(5)*log(T)+chem(c).antConst(6)*T^(chem(c).antConst(7)));
    %chem(c).Psat = @(T) min(chem(c).antConst(1) + chem(c).antConst(2)/(T+chem(c).antConst(3)) + chem(c).antConst(4)*T + chem(c).antConst(5)*log(T)+chem(c).antConst(6)*T^(chem(c).antConst(7)),chem(c).antConst(8);
    chem(c).dHvap = @(T) (chem(c).dHvCon(1)*((1-T/chem(c).Tc)/(1-chem(c).dHvCon(2)/chem(c).Tc))^(chem(c).dHvCon(3)+chem(c).dHvCon(4)*(1-T/chem(c).Tc))/1);
    
    %allow for different inputs for Cp constants
    if isfield(chem(c),'CPGConst')
        if length(chem(c).CPGConst)>5
            chem(c).CPG = @(T) (T>chem(c).CPGConst(7))*(chem(c).CPGConst(1) + chem(c).CPGConst(2)*(T) + chem(c).CPGConst(3)*T^2 + chem(c).CPGConst(4)*(T)^3 + chem(c).CPGConst(5)*T^(4) + chem(c).CPGConst(6)*T^(5))...
                +(T<chem(c).CPGConst(7))*(chem(c).CPGConst(9)+chem(c).CPGConst(10)*T^chem(c).CPGConst(11));
        elseif length(chem(c).CPGConst)>2
            chem(c).CPG = @(T)(chem(c).CPGConst(1) + chem(c).CPGConst(2)*(T) + chem(c).CPGConst(3)*T^2 + chem(c).CPGConst(4)*(T)^(-2));
        else
            chem(c).CPG = @(T) (chem(c).CPGConst(1));
        end
    else
        temp =length(sscanf(chem(c).Formula,'%c %*[0,1,2,3,4,5,6,7,8,9]'));
        chem(c).CPG = @(T)(((5+(temp>1))/2) * 8.3145);
    end
    
    if isfield(chem(c),'CPLConst')
        if length(chem(c).CPGConst)>5
            chem(c).CPL = @(T) chem(c).CPLConst(1) + chem(c).CPLConst(2)*(T) + chem(c).CPLConst(3)*T^2 + chem(c).CPLConst(4)*(T)^3 + chem(c).CPLConst(5)*T^(4);
        elseif length(chem(c).CPGConst)>2
            chem(c).CPL = @(T) chem(c).CPLConst(1) + chem(c).CPLConst(2)*(T) + chem(c).CPLConst(3)*T^2 + chem(c).CPLConst(4)*(T)^(-2);
        else
            chem(c).CPL = @(T) (chem(c).CPLConst(1));
        end
    else
        chem(c).CPL = @(T) (1.2*chem(c).CPG(T));    
    end
    
%     chem(c).CPMeanH = @(T1,T2) integral(chem(c).CPG,T1,T2)/(T1-T2);
%     chem(c).CPMeanS = @(T1,T2) integral((chem(c).CPG/T),T1,T2)/(log(T2/T1));


    %parse SDS
    %%%
%Parse the Hcodes from the SDS
%%%
    
    if isfield(chem(c), 'Hcode')
        hflag=~isempty(chem(c).Hcode);
    else
        hflag=0;
    end
    
    chem(c).hflag=hflag;
    
    if hflag
        hcodes = sort(chem(c).Hcode);
        chem(c).hval=SDSparse(hcodes);
    end

            
end

% Go through CRW report

fid = fopen(crw);
CW= textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

CRAW=deblank(CW{1});
CRAW=strtrim(CRAW);

doci = findLine(CRAW,'DOCUMENTATION REPORT');
hazi = findLine(CRAW,'HAZARDS REPORT');

%split up the documentation and hazard reports. 
cin{1}=1;
cin{2}=1;
in=1;
while in<length(CRAW)
    [B I]=sort([hazi,doci,in]);
    type=find(I==3,1)-1;
    CRW{type}{cin{type},1}=CRAW{in};
    cin{type}=cin{type}+1;
    in=in+1;
end
chem(1).CRW=CRW;

%find which chemicals are in the mixture and put them in 'Matlab' order
chi=findLine(CRW{1},'Chemicals and Reactive Groups in this Mixture')+1;

che=findLine(CRW{1},'------')-1;
che=che(end);
matNames={chem.altName};
for ci= chi:che
    realind= strcmp(matNames,CRW{1}{ci});
    CRWNames{realind,1}=CRW{1}{ci};
end

%go through the hazard report one section at a time.

pointWords{1}={'corrosive','flammable','toxic'};
pointWords{2}={'Intense','pressurization','violent'};
pointWords{3}={'particularly'};
pointVal=[.1,.15,.25,0];
brks= findLine(CRW{1},'--- ');
bgas= findLine(CRW{1},'GASES:');

hin = che+3;
while hin<length(CRW{1})-2
    last=brks(find(brks>hin,1));
    stgas=bgas(find(bgas>hin,1));
    string=strrep(CRW{1}{hin},'mixed with','*');
    string=strrep(string,' -','');
    parse=textscan(string,'%s','delimiter','*');
    parse=parse{1};
    parse = strtrim(parse);
    
    
    flagMult = isempty(strfind(parse{2},'itself'));
    
    %figure out what chemicals are paired together
    
    for go= 1:flagMult+1
        pair(go)=find(strcmp(CRWNames,parse{go}));
    end
    
    temp=hin+2;
    treac=0;
    tgas=0;
    if flagMult
        pair=sort(pair);
    else
        pair=[pair,pair]; 
    end
    
    for hin=temp:last-1
        if ~isempty(strfind(CRW{1}{hin},'No '))
            if hin<stgas
                treac=0;
            else
                tgas=0;
            end
            
            
        else
            if hin<stgas
                %do reactivity hazards
                parse=textscan(CRW{1}{hin},'%s','delimiter',' ');
                parse=parse{1};
                %check thorugh the point values
                got=4;
                for pi=1:length(parse)
                    for gp=1:3
                        if ~isempty(findLine(pointWords{gp},parse{pi}))
                            got =gp;
                        end
                    end
                end
                
                
                    treac=treac+pointVal(got);
                    
                
            else
                if tgas<.5 && hin~=stgas
                
                    tgas=tgas+.1;
                end
                
            end
       
        end
    end
    
    
    ReacMat(pair(1),pair(2))=treac;
    GasMat(pair(1),pair(2))=tgas;

    hin=hin+3;
    
end

chem(1).ReacMat=ReacMat;
chem(1).GasMat=GasMat;





end

% function [chem, exclude] = readHoriz(varargin)
% exclude=[];
% 
% for items=2:length(varargin)
%     ind=locField(fields,varargin{items});
%     exclude=[exclude,ind];
%     cons=1;
%     while ~isnan(raw{ind,indVal-1+cons})
%         chem(c).(fields{ind})(cons,1)=raw{ind+2,indVal-1+cons};
%         cons=cons+1;
%     end
% end
% end

function [ind] = locField(f,string)
ind=find(cellfun(@isempty,strfind(f,string))==0,1);
end

function [out] = locAtom(name,atom)
ind=strfind(name,atom);

if ~isempty(ind)
    result=textscan(name(ind:end),[atom ' %f %*s']);
    if isempty(result{1})
        out=1;
    else
        out=result{1};
    end
else
    out=0;
end

end

function [hval] = SDSparse(hcodes)

hin=1;
    hc = hcodes(1);
    hval=zeros(13,1);

%Safety
    if hc<206
        %Explosion
        hval(6)= abs(hc-206)/6;
        [hc, hin]=update(hcodes,hin);
    end
    if hc< 240
        %Flamability 2
        if hc==226 
            cat=3;
        elseif hc==227
            cat=4;
        elseif (hc/2)==ceil(hc/2)
            cat=1;
        else
            cat=2;
        end

        hval(2)=1.25-cat/4;

        [hc, hin]=update(hcodes,hin);

    end
    if hc<250
        %Expected damage 5
        hval(5)=abs(hc-243)/3;
        [hc, hin]=update(hcodes,hin);

    end
    if hc==250
        %Pyrophoric liquids or solids
        hval(5)=1;
        [hc, hin]=update(hcodes,hin);
    end
    if hc<260
        %self heating
        hval(4)=abs(hc-253)/2;
        [hc, hin]=update(hcodes,hin);
    end
    if hc<270
        %emit flamable gases w/H20
        flagNoH20=1;
        if hc==262
            hval(2)=.9;
        else
            hval(2)=.25;
        end
        [hc, hin]=update(hcodes,hin);
    end
    if hc<280
        %Oxidizer
        hval(2)=abs(hc-273)/3;
        [hc, hin]=update(hcodes,hin);
    end
    if hc<290
        %Gas under pressure
        hval(1)=1;
        [hc, hin]=update(hcodes,hin);
    end
    if hc==290
        %corrosive
        flagCor=1;
        [hc, hin]=update(hcodes,hin);
 %TOX       
    end
    if hc<304
        %acute tox, Oral
        hval(1)=.25;
        hval(3)=abs(hc-304)/5;
        [hc, hin]=update(hcodes,hin);

    end
    if hc<310
        %acute tox swallowed
        hval(3)=.6+abs(hc-306)/5;
        hval(10)=.4+abs(hc-306)/5;
        [hc, hin]=update(hcodes,hin);

    end
    if hc<314
        %acute tox, dermal

        hval(3)=abs(hc-314)/5;
        if hc==310
            hval(3)=1;
        end
        [hc, hin]=update(hcodes,hin);

    end
    if hc<318
        %Irritation
        hval(7)=abs(hc-318)/5;
        if hc==314
            hval(7)=1;
        end
        [hc, hin]=update(hcodes,hin);

    end
    if hc<321
        %Eye Irritation
        hval(7)=abs(hc-321)/5;
        if hc==318
            hval(7)= 1;
        end
        [hc, hin]=update(hcodes,hin);

    end
    if hc<334
        %acute tox Inhaled
        hval(3)=.5+abs(hc-334)/8;
        hval(10)=.4+abs(hc-334)/8;
        [hc, hin]=update(hcodes,hin);



    end
    if hc<340
        %chronic senistization 
        hval(7)=abs(hc-340)/6;
        if hc==334
            hval(7)= 0.8;
        end
        [hc, hin]=update(hcodes,hin);

    end
    if hc<370
        %Carconogins
        hval(8)=.75*(hc/2==ceil(hc/2))+.25;
        [hc, hin]=update(hcodes,hin);

    end
    if hc<400
        %target organs
        cat=5*(hc>371)+3;
        hval(cat)=.65*(hc/2==ceil(hc/2))+.35;
        [hc, hin]=update(hcodes,hin);

 %ENV
    end
    if hc< 410
        %Aquatic, acute
        hval(9)=abs(hc-403)/4;
        [hc, hin]=update(hcodes,hin);
    end
    if hc<420
        %Aquatic, chronic
        hval(9)=abs(hc-414)/5;
        hval(12)=abs(hc-414)/6;
        [hc, hin]=update(hcodes,hin);
    end
    if hc==420
        %Ozone destroyer

        hval(13)=.7;
        hval(10)=.8;
    end



end

function [val, ind] = update(nums , ind)
%new number for the Rcode parsing sequence
if (ind+1) > length(nums)
    val = 999;
    
else
    ind=ind+1;
    val= nums(ind);
end


end