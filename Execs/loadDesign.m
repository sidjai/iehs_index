function [u , s ] = loadDesign(cfg )
%loadDesign Parses the Aspen output to get process conditions and unit info
%   u = unit info in struc array format
%   s = stream info in struc array format
%   p = summary of overall plant connectivity

%xlsname = [cfg.AspenStreams];
%cfg.AspenBlocksName = [cfg.AspenBlocksName];
% clear all;


totdivMean = {'T','P','x','Phase','loc','time'};
divMean = {{'T','rho','CP'},{'P'},{'xd','Activity','rho'},'','','',''};
possChg = {'T','rho','CP','P','xd','Activity'};
startLine = [2;4;2;2;3;2;3];


% [junk junk rawMat] = xlsread(cfg.AspenStreams,1);

findLine = @(c , string) (find(~cellfun('isempty', strfind(c,string))));
%the variables that we want to extract from the xcel sheet Stream
careVars={'Temperature','Pressure','Frac','Flow','Enthalpy', 'Density', 'MW','CP','GAMMA','FIGMX','TBUB','TDEW','DGMX','MUMX'};
fields={'T','P','x','F','H','rho','MW','CP','Activity','fugacity','Tbub','Tdew','dGmix','viscous','Phase'};
% labels{:,1}=txt{:,1};
getVind = @(str) (findLine(careVars,str));

%----------------------------------------------------------------
%----------------------------------------------------------------
%Stream table
%----------------------------------------------------------------
%----------------------------------------------------------------\

rawStm = xlsReadPretty(cfg.AspenStreams,1,'trim');
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
        
        careRef(fr) = length(fields);
    end
end

phaseSet = find(careRef==length(fields))';

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
    pCnt = zeros(length(fields),1);
    
    for asp = allSet(filledSet(:,w))
        id = careRef(asp);
        ref = fields{id};
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
    
    
    
    
end
    
%----------------------------------------------------------------
%----------------------------------------------------------------
%assign the unit parameters (BLOCKS) TXT
%----------------------------------------------------------------
%----------------------------------------------------------------



fid = fopen(cfg.AspenBlocksName);
C= textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

report=strtrim(C{1});

last=length(report);

%open the class designation xcel sheet that assigns the models to classes
rawDes = xlsReadPretty(cfg.ClassDesName,1);
strow = findLine(rawDes{1}(1:end,1),'Model')+1;
vesSet = strow:size(rawDes{1},1);

vesInfo = struct;

for vs = vesSet
    dat = cell2mat(rawDes{1}(vs,[3,5:end]));
    modelType = rawDes{1}{vs,1};
    vesInfo.(modelType) = {dat(1),dat(2),dat(3:end)};
end
    
% [junk models junk] = xlsread(cfg.ClassDesName,1,'A1:A50');
% [cls] = xlsread(cfg.ClassDesName,1);
% [failMat]= xlsread(cfg.ClassDesName,2);

%define the boudaries for the units and parse identifying info
begi=findLine(report,'BLOCK');

%get rid of the instances of errors (has blocks in the line)
beg = delLine(report,begi,'ERRORS',0);
beg = delLine(report,beg,'IMBALANCE',0);

k=1;
for ww=beg'
    parse = textscan(report{ww},'BLOCK: %s MODEL: %s');
    
    u(k,1).name=parse{1}{1};
    u(k,1).type=parse{2}{1};
    
    %asign models
    infoCell = vesInfo.(u(k,1).type);
    u(k,1).failbase = infoCell{1};
    u(k,1).classFunc = infoCell{2};
    u(k,1).classVar = infoCell{3}(~isnan(infoCell{3}));
    
    
    %get the report to the next line
    if k == length(beg)
        fin = last;
        
    else
        fin = beg(k+1)-1;
    end
    ie=1;
    for h=ww:fin
        u(k,1).report{ie,1}=report{h};
        ie=ie+1;
    end
    k=k+1;
end

untLen=k-1;

for n=1:untLen
    
    % first make a list of all the variables that are changing due to the
    % process unit type ie its class.
%     divMean = {{'T','H','rho'},'P',{'xd','Activity','rho'}'','Phase','','',''};
    
    
    
    %get the flowrate of the vessel
    stuff = findLine(u(n).report,'TOTAL BALANCE')+1;
    
    pmole = textscan(u(n).report{stuff},'MOLE (%s)%f%f%f%f');
    pmass = textscan(u(n).report{stuff+1},'MASS (%s)%f');
    penth = textscan(u(n).report{stuff+2},'ENTHALPY (%s)%f%f');
    
    if ~isempty(pmole{5})
        u(n).dF= pmole{4};
        u(n).Fout = [pmole{3}, pmass{2}];
    end
    
    
    %get generation and ouput terms if its a reactor
%     if u(n).classFunc==4
%         pars = textscan(u(n).report{stuff},'MOLE (%s)%f%f%f');
%         u(n).dF = pars{4};
%         u(n).Fout= [pars{3}, pmass{2}];
%     else
%         pars = textscan(u(n).report{stuff},'MOLE (%s)%f');
%     end
    u(n).Fn=pmole{2};
    u(n).Fm=pmass{2};
    
    u(n).F = [pmole{2}; pmass{2}];
%     u(n).report{stuff+1}
%     pmole{1}
%     pmass{1}
    u(n).unitF = {pmole{1}{1} ; pmass{1}{1}};
    
    u(n).dH = penth{3}-penth{2};
    u(n).unitdH = penth{1}{1};
    
    % get the stream I/O info
    in = findLine(u(n).report,'INLET');
    out = findLine(u(n).report,'OUTLET');
    
    in = delLine(u(n).report,in,':',1); %has to have a colon
    in = delLine(u(n).report,in,'PRES',0); %can't be a pressure
%     in = delLine(u(n).report,in,'TEMP',0);%can't be a temp
    
    out = delLine(u(n).report,out,':',1); %has to have a colon
    out = delLine(u(n).report,out,'PRES',0); %can't be a press
%     out = delLine(u(n).report,out,'TEMP',0);%can't be a temp
    
    
    
    lines = [in out'];
    colon = strfind(u(n).report{in},':');
    if length(out)==1
        clout = strfind(u(n).report{out},':');
        things = [colon clout];
        flagSep=0;
    else
        numtop = strfind(u(n).report{out(1)},':');
        numbot = strfind(u(n).report{out(2)},':');
        things = [colon numtop numbot];
        flagSep=1;
    end
    %match up the stream label with the stream number
    
    nin=1;
    flagMix=0;
    for th = things
%         pare{nin} = textscan(u(n).report{lines(nin)}(th+1:end),'%s');
        rem=u(n).report{lines(nin)}(th+1:end);
        li=[];
        while ~isempty(rem)
            [tok, rem] = strtok(rem);
%             li = [li findLine(labels, pare{nin}{1}{pr})];
            
%             findLine(labels, tok)
            li = [li find(strcmpi(labels, tok),1)];
        end
            
        
%         for pr=1:length(pare{nin}{1})
%             u(n).report{lines(nin)}
%             
%         end

        if nin==1
            u(n).in= li; %set the units plant structure
            if length(li)>1
                flagMix=1;
            end
            for er=li
                if estVar(s(er),'to')
                    s(er).to=n; %set the streams too
                end
            end
        else %for the out variables
            for er=li
                if estVar(s(er),'from')
                    s(er).from=n;
                end
            end
            if flagSep
                
                if nin==2
                    u(n).outTop = li;
                else
                    u(n).outBot = li;
                end
            
                
            elseif length(li)>1
                flagSep=1;
                if s(li(1)).Phase>s(li(2)).Phase
                    u(n).outTop = li(1);
                    u(n).outBot = li(2);
                elseif s(li(1)).rho(1)<s(li(2)).rho(1)
                    u(n).outTop = li(1);
                    u(n).outBot = li(2);
                else
                    u(n).outTop = li(2);
                    u(n).outBot = li(1);
                end
                
            else
                u(n).out = li;
                
            end
        end
        nin=nin+1;
        
    end
    u(n).multOut = flagSep;
    u(n).multIn = flagMix;
    
    if flagSep
        u(n).out =[u(n).outTop, u(n).outBot];
    end
    
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
            if ~isempty(divMean{temp(k)})
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
                chu=divMean{d};
                chk=[];
                for fc=1:length(chu)
                    chk(fc)=estVar(u(n),chu{fc});
                end
                if max(chk)==1
                    des{1,k} = [des{1,1} ' w/ ' totdivMean{d} '= ' unitChg{a}];
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
for su=1:length(s)
    if estVar(s(su),'from')
        s(su).from=0;
    elseif estVar(s(su),'to')
        s(su).to=999;
    end
end


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

function [ s ] = loadAspenStreams(fullfile)




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
