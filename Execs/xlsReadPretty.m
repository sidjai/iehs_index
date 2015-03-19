function [ out ] = xlsReadPretty(varargin)
%xlsReadPretty Read xls documents by converting to csv
%   numSheetsParsed is the second var, if it is not provided the program
%   the next args are options for postprocessing
%   if windows, uses vbscript
%   if linux/mac, searches for libreoffice using unoconv --format  xls
%   example.csv
fullfile = varargin{1};
argLen = length(varargin);
if argLen == 1
    numSheetsParsed = 9999;
else
    numSheetsParsed = varargin{2};
end

if argLen <= 2 
    stdin = {};
else
    stdin = varargin(3:end);
end

realWd = strrep(pwd,'\','/');
comd = [realWd '/xls2csv.vbs "' fullfile '" ' int2str(numSheetsParsed)];
[ranTest] = system(comd);
if(logical(ranTest))
    error('vbscript messed up')
end

ind = regexp(fullfile,'[.]');
beg = fullfile(1:(ind(end)-1));

sh = 1;
csvFile = [beg num2str(sh) '.csv'];
while(exist(csvFile, 'file'))
    fid = fopen(csvFile,'r');
    parse = textscan(fid, '%s','delimiter','\n');
    parse = parse{1};
    fclose(fid);
    for (row = 1:length(parse))
        line = [parse{row} ','];
        commas = regexp(line,',');
        col = 2;
        bef = commas(1);
        if bef~=1
            out{sh}{row,1} = line(1:bef-1);
        end
        
        for ca = commas(2:end)
            %if the commas are next to each other just add col
            
            if (bef+1 ~= ca)
               out{sh}{row,col} = line(bef+1:ca-1);
               
            end
            col=col+1;
            bef = ca;
        end
        
    end
    
    %clean up
    %need numbers not strings so convert those that need it
    out{sh} = doOption(out{sh},'EmptyisNaN');
    out{sh} = doOption(out{sh},'OneSpaceIsEmpty');
    out{sh} = doOption(out{sh},'NumsNotStrings');
    %Do the optional utilities
    for op = 1: length(stdin)
        out{sh} = doOption(out{sh},stdin{op});
    end
    delete(csvFile)

    sh = sh+1;
    csvFile = [beg num2str(sh) '.csv'];
    
end


end
function [cin] = doOption(cin,optionTag)

switch optionTag
    case 'trim'
        badSet = cellfun(@(x)ischar(x),cin);
        cin(badSet) = strtrim(cin(badSet));
        
    case 'NumsNotStrings'
        badSet = cellfun(@(x)...
           ischar(x) && ~isempty(x)...
           && isempty(regexpi(x,'[a-df-z]','once'))...
           ,cin);
        cin(badSet) = cellfun(@(x)str2double(x),...
            cin(badSet),'UniformOutput',false);
    case 'EmptyisNaN'
        badSet = cellfun(@(x)isempty(x),cin);
        cin(badSet) = {NaN};
    case 'OneSpaceIsEmpty'
        oneSet = cellfun(@(x)( length(x) == 1 && ~isnan(x)),cin);
        cin(oneSet) = strrep(cin(oneSet),' ','');
end


end
