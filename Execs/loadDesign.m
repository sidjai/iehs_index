function [u , s ] = loadDesign(cfg )
%loadDesign Parses the Aspen output to get process conditions and unit info
%   u = unit info in struc array format
%   s = stream info in struc array format
%   p = summary of overall plant connectivity

%xlsname = [cfg.AspenStreams];
%cfg.AspenBlocksName = [cfg.AspenBlocksName];
% clear all;


% totdivMean = {'T','P','x','Phase','loc','time'};
% divMean = {{'T','rho','CP'},{'P'},{'x','Activity','rho'},'','','',''};
% possChg = {'T','rho','CP','P','x','Activity'};
% startLine = [2;4;2;2;3;2;3];


% [junk junk rawMat] = xlsread(cfg.AspenStreams,1);

findLine = @(c , string) (find(~cellfun('isempty', strfind(c,string))));
%the variables that we want to extract from the xcel sheet Stream

vesselDict = getClassDes(cfg.ClassDesName);
classVarDef = vesselDict.classVarDef;
allchg = classVarDef(:,2);

k=1;
for rw = 1:length(allchg)
    cChg = allchg{rw};
    for co = 1:length(cChg)
        possChg{k,1} = cChg{co};
        k=k+1;
    end
end

possChg = unique(possChg); %stop alphabetical order
s = loadStreams(cfg.StreamTable);

if (cfg.conceptDesFlag)
    u = loadConceptBlocks(cfg.ConceptBlock,vesselDict,s);
    
else
    u = loadAspenBlocks(cfg.AspenBlocksName,vesselDict ,{s.name});
end

untLen = length(u);

%Do the overarching structure attributions, streams attaching to what
transIn = {u.in};
transOut = {u.out};
mixSet = cellfun(@(x)length(x)>1,transIn);
splitSet = cellfun(@(x)length(x)>1,transOut);



for ms = 1:untLen
    u(ms).multIn = mixSet(ms);
    u(ms).multOut = splitSet(ms);
    for sind = transIn{ms}
        s(sind).to = ms;
    end
    
    for sind = transOut{ms}
        s(sind).from = ms;
    end
end

prodSet = cellfun(@(x)isempty(x),{s.to});
feedSet = cellfun(@(x)isempty(x),{s.from});
for k = find(feedSet)'
    s(k).from = 0;
end

for k = find(prodSet)'
    s(k).to = 999;
end

for fk = find(splitSet)'
    test = [u(fk).out',...
        [s(u(fk).out(1)).Phase(1);s(u(fk).out(2)).Phase(1)]];
    test = sortrows(test,2);
    u(fk).outTop = test(1,1);
    u(fk).outBot = test(2,1);
end

for n=1:untLen
    
    
   %parse the unit information
   want{2,1}={'LATE TEMP','OM TEMP','NET CON' , 'NET REB','ACTUAL EQUI','HEAT DUTY','TO FEED RATIO'}; %Mass transfer units
   want{3,1}={'HEAT DUTY'};%Physical changes-T 'heater'
   want{3,2}={'LET PRES','LET TEMP','IC EFF','BRAKE','HEAD D','HEAT CAPA'};%Physical changes-P 'compressor'
   want{4,1}={'HEAT DUTY','SYSTEM TEMP','SYSTEM PRES','PHASE FRAC','CONV FRAC','   EXTENT','RY MATRIX','T SPEC'};%Reactor
   want{5,1}={};%transport
   
   namesC{2,1}={'Ttop','Tbot','dutyCon','dutyReb','stages','duty','refluxRatio'};
   namesC{3,1}={'duty'};
   namesC{3,2}={'P','T','effIsen','power','Head','gamma'};
   namesC{4,1}={'duty','T','P','Phase','conv','ci','stoich','inert'};
   namesC{5,1}={};
   
   wind= u(n).classFunc;
   if wind ==3
       u(n).flagCompr = u(n).classVar(1)-1;
       wind=sub2ind(size(want),wind,u(n).classVar(1)); %if the first changed variable is T then use {3,1}
       
       
   end
   
   [u]= parseReport(u, n,want{wind},namesC{wind});
   
%     divMean = {{'T','H','rho','CP'},{'P'},{'xd','Activity','rho'},'','','',''};
% possChg = {'T','H','rho','CP','P','xd','Activity'};

    
    %find out which variable changes
    
    possLen= length(possChg);
    out=u(n).out;
    in=u(n).in;
    chg=zeros(possLen,1);
    flagChgInMix=chg;
    
    %do the base case
    for b=1:possLen
        if estVar(u(n),possChg{b})
            chg(b,1)=out(1);
        end
    end
    unitBase='Outlet';
    
    
    if u(n).classFunc<5
        temp=u(n).classVar;
        j=1;
        for k=1:length(temp)
            if ~isempty(classVarDef{temp(k),2})
                divCh(j)=temp(k);
                j=j+1;
            end
        end
        
        
        %find out how many alternatives you need
        %inlet, oulet for each of the changes, approx avg for everything
        
        numAlt=1+length(out)*length(divCh)+1;
        
        flagChgInMix=zeros(possLen,numAlt+1);
        
        
        
        if u(n).multOut
            unitBase='Avgout';
            sBase =-888; %average the two outlets
            unitChg={'Inlet','Top','Bot'};
            schg =[in(1),out(1),out(2)];
            
        elseif u(n).classFunc == 3
            
            unitBase='Avg';
            sBase = -999; %flag to avg the in and out
            unitChg ={'Inlet','Outlet'};
            schg= [in(1),out(1)];
            
        else
            
            unitBase='Outlet';
            sBase = out;
            unitChg ={'Inlet'};
            schg= in(1);
        end
        
        
        for b=1:possLen
            if estVar(u(n),possChg{b})
                chg(b,1)=sBase;
            end
        end
        allSet = find(chg(:,1))';
        
        des{1,1} =['Base:' unitBase];
        k=2;
        for a=1:length(schg)
            for d=divCh
                chu = classVarDef{d,2};
                chk=[];
                for fc=1:length(chu)
                    chk(fc)=estVar(u(n),chu{fc});
                end
                if max(chk)==1
                    des{1,k} = [des{1,1} ' w/ ' classVarDef{d,1} '= ' unitChg{a}];
                    for ff=1:length(chu)
                        if chk(ff)
                            nu=find(strcmp(possChg,chu{ff}),1);
                            chg(nu,k)= schg(a);
                        end

                        if a==1;
                            flagChgInMix(nu,k)=u(n).multIn;
                        end

                    end
                    chg(chg(allSet,k)==0,k)=sBase;
                    k=k+1;
                end
            end

            for p=find(chg(:,1))
                chg(p,k)=schg(a);
                des{1,k} = ['Base:' unitChg{a}];
                k=k+1;
            end
                
    
            
        end
                
            
     
        
    
        
    else
        des{1,1} =['Base:' unitBase];
    end
        
    
    u(n).chg=chg;
    
    u(n).flagChgInMix=flagChgInMix;
    u(n).maDes=des;
    
    u(n).numCon=size(chg,2);
    u(n).chgFlag=u(n).numCon>1;

    

    
    
    
    
%         u(n).V=(u(n).stages*70) * 
        
 
        
        
   
%         s(u(n).out).(potFields{18})
%         u(n).out
%     end
        
    
    
   
    
    
end

%Go back and do more complete labels
for n=1:untLen
    %Add a name for the input and output on the unit itself
    IO{1}=u(n).in;
    IO{2}=u(n).out;
    IOval{1}='Input: ';
    IOval{2}='Output: ';
    for io=1:2

        for pi=1:length(IO{io})
            if pi<length(IO{io})
                IOval{io}=[IOval{io} s(IO{io}(pi)).name '; '];
            else
                IOval{io}=[IOval{io} s(IO{io}(pi)).name];
            end
        end
    end
    u(n).Indes=IOval{1};
    u(n).Outdes=IOval{2};
    
    
end




%Go back through the streams and assign the streams with feed and product



end

function [val] = findLine(c,st)
val =(find(~cellfun('isempty', strfind(c,st))));
end

function [u] = parseReport(u,n,ca, name)
% for the specified unit, go into the report and find the variables in ca
% then assign the value reported as a field in the structure for the unit.
% Also assigns the units for the values as unit"val" ie. unitT.

for er = 1: length(ca)
    line =findLine(u(n).report,ca{er});

    if ~isempty(line)

        stemp=u(n).report{line(1)};
        stemp=strrep(stemp,',','');

        if ~isempty(strfind(stemp,'('))
            %if it has the units in parenthesis, shortcut

            parse = textscan(stemp,'%s','delimiter','()');
            u(n).(name{er}) = str2num(parse{1}{3});
            u(n).(['unit' name{er}]) = strtrim(parse{1}{2});
        elseif u(n).classFunc==4
            if ~isempty(strfind(stemp,'CONV FRAC'))
                for r=1:length(line)
                    stemp=u(n).report{line(r)};
                    nums = strfind(stemp,':');
                    parse = textscan(stemp,'%s','delimiter',':');
                    u(n).(name{er})(r,1)=str2num(parse{1}{end});
                    u(n).(['unit' name{er}]){r,1} = strtok(parse{1}{3});
                end


            elseif ~isempty(strfind(stemp,'EXTENT '))
                line=line+1;
                u(n).(['unit' name{er}]) = u(n).report{line};
                line=line+1;
                stemp=u(n).report{line};
                while ~isempty(stemp)
                    
                    parse=textscan(stemp,'%f%f');
                    u(n).(name{er})(parse{1},1)=parse{2};
                    line=line+1;
                    stemp=u(n).report{line};
                end

            elseif ~isempty(strfind(stemp,'MATRIX'))
                % get the stoic matrix, 
                %chem name    number   chem name ...
                line=line+2;
                stemp=u(n).report{line};
                hp=strfind(stemp,'#');
                while ~isempty(hp)
                    hparse = textscan(stemp(hp(end):end),'#%f');
                    r=hparse{1};
                    line=line+2;
                    stemp= u(n).report{line};

%                     stoic=zeros(length(chem),1);
                    fl=1;
                    [tok rem]=strtok(stemp);
                    while fl>1;

                        u(n).stoCel{r}{fl,1}=tok;
%                         c=findLine({chem.ID},tok);
                        [tok rem]=strtok(rem);
                        u(n).stoCel{r}{fl,2}=str2num(tok);

                        if isempty(rem)
                            fl=0;
                        else
                            [tok rem]=strtok(rem);
                        end

                    end



                    line=line+2;
                    stemp=u(n).report{line};
                    hp=strfind(stemp,'#');

                end
            elseif ~isempty(strfind(stemp,'SPEC'))
                line =line+1;
                stemp=u(n).report{line};
                fl=1;
                while ~isempty(stemp)
                    parse = textscan(stemp,'%s%f%s');
                    if fl==1;
                        u(n).unitinert=parse{3};
                    end
                    
                    u(n).inert{fl,1}=parse{1};
                    u(n).inert{fl,2}=parse{2};
                    fl=fl+1;
                    line=line+1;
                    stemp=u(n).report{line};
                end
            else
                p=getNums(stemp);
            
                parse = textscan(stemp(p(end):end),'%f');
                u(n).(name{er}) = parse{1};

                if length(p)>1
                     u(n).(['unit' name{er}]) = strtrim(stemp(p(end-1):p(end)));
                end
            end
                
                        


        else

            %Do it the fail safe way by separating by two spaces
%                 parse = textscan(stemp, '%s');
            p=getNums(stemp);
            
            parse = textscan(stemp(p(end):end),'%f');
            u(n).(name{er}) = parse{1};
            
            if length(p)>1
                 u(n).(['unit' name{er}]) = strtrim(stemp(p(end-1):p(end)));
            end

        end
    end
end
end

function [cnum] = getNums(st)
nums = strfind(st,'  ');
if isempty(nums)
    cnum=0;
else
    k=2;
    cnum(1) = nums(1);
    for ni=2:length(nums)
        if ~(1==(nums(ni)-nums(ni-1)))
            cnum(k)=nums(ni);
            k=k+1;
        end
    end
end

end

function [ s ] = loadStreams(fullfile)
%Load the stream list from an output of Aspen


careVars = {'Temperature','Pressure','Frac','Flow','Enthalpy', 'Density', 'MW','CP','GAMMA','FIGMX','TBUB','TDEW','DGMX','MUMX'};

careFields={'T','P','x','F','H','rho','MW','CP','Activity','fugacity','Tbub','Tdew','dGmix','viscous','Phase'};
% labels{:,1}=txt{:,1};
getVind = @(str) (findLine(careVars,str));


rawStm = xlsReadPretty(fullfile,1,'trim');
%load up margins with the text on the left and headings for the top stuff
rawStm = rawStm{1};
margin = rawStm(:,1);
headLen = find(cellfun('length', margin)>1,1);
stmLen = size(rawStm,2)-1;
margin = margin(headLen+1:end);
valTable = rawStm(headLen+1:end,2:end);

%load up the headings for the streams
headings = rawStm(1:headLen,2:stmLen+1);
for col= 2: (stmLen+1)
    s(col-1).name = rawStm{1,col};
    s(col-1).labels  = rawStm(1:headLen,col); 
end
careRef = zeros(length(margin),1);
unitAsp = cell(length(margin),1);
valuAsp = cell(length(margin),stmLen);

%----------------------------------------------------------------
%Figure out what the margins are refering to in terms of fields
%----------------------------------------------------------------
for j = 1:length(margin)
    temp = cellfun(@(x) strfind(margin{j},x),careVars,'UniformOutput',0);
    ind = find(cellfun('length',temp),1);
    if isempty(ind) ...
        || (ind==getVind('Flow') && isempty(strfind(margin{j},'otal')))
        ind = 0;
    end
    careRef(j) = ind;
end
allSet = find(careRef)';
fracSet = find(careRef==getVind('Frac'))';
for fr = fracSet
    if (~isempty(strfind(margin{fr},'apor'))...
        || ~isempty(strfind(margin{fr},'iquid'))...
        || ~isempty(strfind(margin{fr},'olid')))
        
        careRef(fr) = length(careFields);
    end
end

phaseSet = find(careRef==length(careFields))';

vertSet = find((careRef==getVind('Frac') ...
        | careRef==getVind('GAMMA')))';

horzSet = allSet(~ismember(allSet,union(vertSet, phaseSet)));
%Go through different situations and return the value that should be
%returned and the units associated with it.

for c = horzSet
    parsenum = regexp(margin{c},'\s');
    unitAsp{c} =  margin{c}(parsenum(end)+1:end);
end

for dw = union(vertSet, phaseSet)
    unitAsp{dw} = strtok(margin{dw});
end

%units all done
%now do the actual values

%get the values for the entries with multiple vertical entries

for dw = vertSet
    ele = dw+1;
    %get the chemicals by their being no lowercase letters
    while isempty(regexp(margin{ele},'[a-z]', 'once'))
        for w = 1:stmLen
            temp = valTable{ele,w};
            if isempty(temp)
                if ~isempty(strfind(margin{dw},'GAMMA'))
                    temp = 1;
                else
                    temp = 0;
                end
            end
                
            valuAsp{dw,w}(1,ele-dw) = temp;
        end
        ele = ele+1;
    end
    
end

%As we have the number of elements, get the chem names
ci = 1;
for c = dw+1:ele-1
    s(1).chemNames{ci,1} = margin{c};
    ci=ci+1;
end
    

%get the values for the other entries
for rw = union(horzSet, phaseSet)
    for w = 1:stmLen
        valuAsp{rw,w} = valTable{rw,w};
    end
end

%set everything into their places in the structure
filledSet = cellfun(@(x)~isempty(x),valuAsp(allSet,:));
for w = 1:stmLen
    pCnt = zeros(length(careFields),1);
    
    for asp = allSet(filledSet(:,w))
        id = careRef(asp);
        ref = careFields{id};
        pCnt(id) = pCnt(id) +1;
        
        s(w).(ref)(pCnt(id),:) =  valuAsp{asp,w};
        if pCnt(id) == 1
            s(w).(['unit' ref ]) = unitAsp{asp};
        elseif pCnt(id) == 2
            tem = s(w).(['unit' ref ]);
            s(w).(['unit' ref ]) = cell(2,1);
            s(w).(['unit' ref ]){1} =tem;
            s(w).(['unit' ref ]){2} = unitAsp{asp};
        else
            s(w).(['unit' ref ]){pCnt(id)} = unitAsp{asp};
        end
    end
end

%Add some additional info

%Change the CP units to whatever we see fit
cpunits = {s.unitCP};
filled = cpunits(cellfun(@(x)~isempty(x),cpunits));
chSet = findLine(filled,'kg');
for del = filled(chSet)
    s(del).CP=s(del).CP .* s(del).MW(1);
    for ei=1:length(s(del).unitCP)
        s(del).unitCP{ei}=strrep(s(del).unitCP{ei},'kg','kmol');
    end
end

for w = 1:stmLen
    s(w).compList = find(s(w).x(1,:) > 10^-8);
    s(w).numComp = length(s(w).compList);
    
    labels{w,1} = s(w).name;
    
    %make new flow fields that are useful
    s(w).Fn=s(w).F(1);
    s(w).Fm=s(w).F(2);
    s(w).Fv=s(w).F(3);
    
    %Change the 
    
    
    
end



end

function [ u ] = loadAspenBlocks(blockfile, vesInfo, stmNames)
%Grab Aspen blocks text file

fid = fopen(blockfile);
C= textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

fullReport=strtrim(C{1});



%define the boudaries for the units and parse identifying info
dotLine = findLine(fullReport,'-------------------------');

blockEnd = [dotLine(2:end)-2; length(fullReport)];


k = 1;
for lind = dotLine'
    title = cleanCstr(regexp(fullReport{lind-1},'(\w*[:])\s','split'));
    
    u(k,1).name = title{1};
    u(k,1).type = title{2};
    
    %asign models
    infoCell = vesInfo.(u(k,1).type);
    u(k,1).failbase = infoCell{1};
    u(k,1).classFunc = infoCell{2};
    u(k,1).classVar = infoCell{3};
    
    indvRep = fullReport(lind-1:blockEnd(k));
    
    
    %do in out streams in the report with a semicolon
    
    direction = {'in','out'};
    for ty = 1:2
        usedStm = parseStrwMargin(indvRep{ty+2},':');
        u(k,1).(direction{ty}) = cellfun(...
            @(us)find(strcmp(stmNames,us),1)...
            ,usedStm);
    end
    
    %Mass and Energy balances
    baLine = findLine(indvRep,'TOTAL BALANCE')+1;
    domain = {'Fn','Fm','dH'};
    
    balanceMat = cellfun(...
        @(x) parseStrwMargin(x,''),...
        indvRep(baLine:baLine+2),'UniformOutput',0);
    balanceMat{1} = balanceMat{1}(1:4);
    balanceMat = cat(1,balanceMat{:});
    
    %Reaction or change in moles
    if (balanceMat{1,4}>.001)
        u(k,1).dF = balanceMat{1,4};
        u(k,1).Fout = [balanceMat{1:2,3}];
    end
    
    mat = [balanceMat{1:2,2},balanceMat{3,4}];
    
    for ty = 1:3
        u(k,1).(domain{ty}) = mat(ty);
        u(k,1).(['unit' domain{ty}]) = balanceMat{ty,1};
    end
    
    u(k,1).F = [balanceMat{1:2,2}]';
    u(k,1).unitF = balanceMat(1:2,1);
    
    u(k,1).report = indvRep;
    k=k+1;
end




end

function [ dict ] = getClassDes(desName)
%open the class designation xcel sheet that assigns the models to classes
rawDes = xlsReadPretty(desName,2);
strow = findLine(rawDes{1}(1:end,1),'Model')+1;

vesSet = strow:size(rawDes{1},1);

dict = struct;

for vs = vesSet
    dat = cell2mat(rawDes{1}(vs,[3,5:end]));
    modelType = rawDes{1}{vs,1};
    classVar = dat(3:end);
    classVar = classVar(~isnan(classVar));
    dict.(modelType) = {dat(1),dat(2),classVar};
end
strow = find(cellfun(@(x)...
    ischar(x) && strcmp(x,'ClassVar'),...
    rawDes{2}(1:end,1)))+1;
    
%Get the variable change definition from the second page
rowSet = strow:size(rawDes{2},1);
dict.classVarDef = rawDes{2}(rowSet,2);

width = size(rawDes{2},2);
for k = 1:length(rowSet)
    cols = rawDes{2}(rowSet(k),3:width);
    goodSet = cellfun(@(x)ischar(x),cols);
    dict.classVarDef{k,2} = cols(goodSet);
end

end

function [ u ] = loadConceptBlocks(blockfile, vesInfo,s)

rawBlock =  xlsReadPretty(blockfile,1);
rawBlock = rawBlock{1};
strow = findLine(rawBlock(1:end,1),'Blocks')+1;

fieldCol = findLine(rawBlock(1,1:end),'Field');


u(1).name = rawBlock{strow,2};
u(1).type = 'CONCEPT';
extVar =struct('T',[],'P',[]);
rows = strow:size(rawBlock,1);
vars = rawBlock(rows,fieldCol);
res = rawBlock(rows,2);
extSet = cellfun(@(x)isfield(extVar,x),vars);
for ex = find(extSet)'
    res{ex} = res{ex} + [s(1).(vars{ex})];
end

for v = 1:length(vars)
    u(1).(vars{v}) = res{v};
end


u(1).in = 1;
u(1).out = 2;
% 
% u(1).T = s(1).T + rawBlock{strow+1,2};
% u(1).P = s(1).P + rawBlock{strow+2,2};
% u(1).dHreac = rawBlock{strow+3,2};
% u(1).V = rawBlock{strow+4,2};

%asign models
infoCell = vesInfo.(u(1).type);
u(1).failbase = infoCell{1};
u(1).classFunc = infoCell{2};
u(1).classVar = infoCell{3};

end
function [cstr] = cleanCstr ( cin )
cin = cellfun(@(x) strtrim(x),cin,'UniformOutput', 0);
goodSet = cellfun(@(x) ~isempty(x),cin);
cstr = cin(goodSet);
end

function [parse] = parseStrwMargin(str, marg)

%margin string ignored if it has a parenthesis

parse = {};

paren = strfind(str,')');
if ~isempty(paren)
    parse = regexp(str(1:paren),'\((\S*)','tokens');
    parse = parse{1};
    tok = str(paren+1:end);
else
    tok = str(regexp(str,marg,'end')+1:end);
end

add = regexp(tok,'(\S*)','tokens');

parse = [parse,add{:}];
parse = cleanCstr(parse);

%find numbers
numSet = cellfun(@(x)...
           isempty([regexpi(x,'[a-df-z]?') regexpi(x,'[a-df-z0-9]-\S?')])...
           ,parse);
for n = find(numSet)
    parse{n} = str2double(parse(n));
end

end

function [val] = delLine(pot,lines,st,flag)

exclude=findLine({pot{lines}},st);
if flag
    loge=zeros(size(lines));
else
    loge=ones(size(lines));
end
loge(exclude)=flag;
loge=logical(loge);
val=(lines(loge));



end
