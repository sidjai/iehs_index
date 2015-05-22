function [ dat ] = xlsReadPretty(varargin)
%xlsReadPretty Read xls documents by converting to csv
%   numSheetsParsed is the second var, if it is not provided the program
%   the next args are options for postprocessing
%   if windows, uses vbscript
%   if linux/mac, searches for libreoffice using unoconv --format  xls
%   example.csv
fullpath = varargin{1};
argLen = length(varargin);
if argLen == 1
    numSheetsParsed = 9999;
else
    numSheetsParsed = varargin{2};
end

if argLen <= 2 
    stdin = {};
    head = 0;
else
    stdin = varargin(3:end);
    if(isempty(findLine(stdin, 'header')))
      head = 0;
    else 
      head = 2;
    end
end

realWd = strrep(pwd,'\','/');

ind = regexp(fullpath,'[.]');
beg = fullpath(1:(ind(end)-1));


if (isunix())
	comd = ['ssconvert -S "' fullpath '" "' beg '%n.csv"'];
else
	comd = [realWd '/xls2csv.vbs "' fullpath '" ' int2str(numSheetsParsed)];
end

[ranTest] = system(comd);
if(logical(ranTest))
	error('convert xlsx to csv messed up')
end


csvFile = [beg num2str(ispc()) '.csv'];
sh = 1;
while(exist(csvFile, 'file'))
  dat{sh}= csvReadPretty(csvFile, head);
	delete(csvFile)
	
  %clean up
  %need numbers not strings so convert those that need it
  dat{sh} = doOption(dat{sh},'EmptyisNaN');
  dat{sh} = doOption(dat{sh},'OneSpaceIsEmpty');
  dat{sh} = doOption(dat{sh},'NumsNotStrings');
  %Do the optional utilities
  for op = 1: length(stdin)
      dat{sh} = doOption(dat{sh},stdin{op});
  end
  
  sh = sh+1;
	csvFile = [beg num2str(sh - isunix()) '.csv'];
end


end

function [out] = csvReadPretty(csvPath, hd)
	fid = fopen(csvPath,'r');
    parse = textscan(fid, '%s','delimiter','\n');
    parse = parse{1};
    fclose(fid);
    for (row = (hd+1):length(parse))
        line = [parse{row} ','];
        commas = regexp(line,',');
        col = 2;
        bef = commas(1);
        if bef~=1
            out{row-hd,1} = line(1:bef-1);
        end
        
        for ca = commas(2:end)
            %if the commas are next to each other just add col
            
            if (bef+1 ~= ca)
               out{row,col} = line(bef+1:ca-1);
               
            end
            col=col+1;
            bef = ca;
        end
        
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
