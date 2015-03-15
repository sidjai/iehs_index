function [ confirm ] = xlsWritePretty( sheets, names, fullfile)
%xlsWritePretty Write the excel sheets by going through csv
%   Sheets have sheets as seperate cells and are cells with the grid
%   names is a cell with all the names of the sheets


ind = regexp(fullfile,'[.]');
beg = fullfile(1:(ind(end)-1));
csvFile = [beg '.csv'];

realWd = strrep(pwd,'\','/');

for sh = 1: length(sheets)
    
    dms = size(sheets{sh});
    
    
    %make the csv
    convertCell = cellfun(@(x)~(ischar(x) || isempty(x)),sheets{sh});
    sheets{sh}(convertCell) = cellfun(@(x)num2str(x),...
        sheets{sh}(convertCell),'UniformOutput',false);
    %sheets{sh} = cellfun(@(x)[x ','],sheets{sh},'UniformOutput',false);
    
    for row =1:dms(1)
        csvLines{row} = sprintf('%s,',sheets{sh}{row,:});
        %csvLines{row} = [sheets{sh}{row,:}];
    end
    
    %write file
    fid = fopen(csvFile,'w');
    cellfun(@(x) fprintf(fid,'%s\n',x),csvLines);
    fclose(fid);
    
    
    
    %call VBscript
    
    comd = [realWd '/csv2xls.vbs "' fullfile '" ' int2str(sh) ' ' names{sh}];
    [ranTest(sh)] = system(comd);
    if(logical(ranTest(sh)))
        error(['Error: command ' comd ' failed'])
    end
    
    csvLines = [];
    
    
end
    
confirm = max(ranTest);

end

