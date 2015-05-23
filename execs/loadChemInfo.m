function [ chem ] = loadChemInfo( aspenPred, TESTDir )
%loadChemInfo Loads basic information about the chemicals involved in the plant
%design
%   Goes through the Aspen chemical prediction report and the CRW
%   report for useful information

estNames = {'TC','PC','OMEGA','CPIG','PLXANT','DHVLWT','DHFORM','DGFORM'};
estFields = {'Tc','Pc','acentric','CPGConst','antConst','dHvCon','dHform','dGform'};


report = readTextSimple(aspenPred,90);

beg = findLine(report,'COMPONENT ID');
k=1;

%go through the report to get the process conditions
for ww=beg'
    parse = textscan(report{ww},'%s');
    chem(k).ID=parse{1}{3};
    temp = strfind(parse{1}{5},'-');
    if isempty(temp)
        chem(k).Formula=parse{1}{5};
    else
        chem(k).Formula=parse{1}{5}(1:temp-1);
    end



    %get the report to the next line and put it in the specific chem
    if k == length(beg)
        fin = length(report);
    else
        fin = beg(k+1)-2;
    end

    ie=1;
    for h=ww:fin
        chem(k).report{ie,1}=report{h};
        ie=ie+1;
    end
    k=k+1;
end
chmLen = k-1;

%find the estimated parameters that are useful
for c=1:chmLen
    %find where the parameters show up
    iue=1;
    for par = 1:length(estNames)
        temp = findLine(chem(c).report,estNames{par});
        neg = findLine(chem(c).report, 'AT');
        temp = temp(~ismember(temp,neg));
        if ~isempty(temp)
            uwe=strfind(chem(c).report{temp},estNames{par})+length(estNames{par});

            inds(iue,:) =[par, temp, uwe];
            iue=iue+1;
        end
    end

    for in=1:size(inds,1)
        var=estFields{inds(in,1)};

        text=chem(c).report{inds(in,2)};
        ad=1;
        rightc=strfind(text,')');

        %if the parameter is temperature dependent, get all the constants
        while ~isempty(rightc)
            if ad ==1;
                text=text(rightc:end);
                parse = textscan(text,') %f %s %*s');
                chem(c).(['unit' var]) = parse{2}{1};
            else
                parse = textscan(text,') %f %*s');
            end

            chem(c).(var)(ad,1) = parse{1};

            text = chem(c).report{inds(in,2)+ad};
            rightc=strfind(text,')');
            text=text(rightc:end);
            ad=ad+1;
        end
        %if its a pure component just get the variable and units
        if ad==1
            text=text(inds(in,3):end); %get text from the var name to the end
            parse = textscan(text,'%f %s %*s');
            val= parse{1};

            unitName= parse{2}{1};

            if strcmp(unitName,'J/KMOL')
                chem(c).(var) =val/1000;
                chem(c).(['unit' var]) =['K' unitName];
            else
                chem(c).(var) =val;
                chem(c).(['unit' var]) =unitName;
            end
        end

    end


end

%run TEST automatically

%Extract the QSAR properties (TEST software required)

%get all files in the folder with Consens in the title
cont = dir(TESTDir);
rawNm = {cont.name};
consLog = findLine(rawNm,'Consen');
txtLog = findLine(rawNm,'txt');

testNm = { cont(intersect(consLog,txtLog)).name};
testLen = length(testNm);
%figure out what fields each txt file refers to
fieldNm = cell(testLen,1);
%tags in the txt file
matchNm = {'Bioaccum','Daphnia','minnow','Flash','boiling','rat','pyrifor'};
%coresponding field it refers to
fieldsTot = {'BCFpred','MGLC50','FMLC50','Flash','Tb','LCInhal','IGC50'};
%transfers the chemID to the chem number already in the system
% chemDict = [2,3,4,5,7,8];
chemDict = [9,7,12,10,5,1,3,8,4,11,6,2];
for ti = 1:testLen
    fieldNm{ti} = fieldsTot{findMatch(testNm{ti},matchNm)};
end

for ti = 1:testLen
    fid = fopen([TESTDir testNm{ti}]);
    raw = textscan(fid, '%s %s %s %s %s %s',...
        'Delimiter', '\t','HeaderLines',1,'CollectOutput' ,true );
    fclose(fid);
    raw = raw{1};
    ctLen = size(raw,1);
    post = zeros(ctLen,2);
    unit = cell(ctLen);
    %get the literature value or the consensus value post(:,2)
    nolitSet = findLine(raw(:,3),'N/A');
    nopredSet =  findLine(raw(:,4),'N/A');
    allSet = 1:ctLen;
    litSet = allSet(~ismember(allSet,nolitSet));
    predSet = allSet(~ismember(allSet,nopredSet));
    predSet = predSet(~ismember(predSet,litSet));

    %Get the values
    post(predSet,2) = str2double(raw(predSet,4));
    post(litSet,2) = str2double(raw(litSet,3));
    %get the chem numbers
    post(allSet,1) = chemDict(allSet)';

    %Set the units
    for ps = predSet
        unit{ps} = 'Predicted';
    end
    for lis = litSet
        unit{lis} = 'Literature';
    end


    %set the cell table in the right field
    for ci = find(post(:,2))'
        chem(post(ci,1)).(fieldNm{ti}) = post(ci,2);
        chem(post(ci,1)).(['unit' fieldNm{ti}]) = unit{ci};
    end


end

end

function [val] = findMatch(str,cell)
val=0;
k=1;
while val == 0 && k <=length(cell)
    temp = strfind(str,cell{k});
    if isempty(temp)
        k = k+1;
    else
        val = k;
    end

end
end
