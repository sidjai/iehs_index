function [ confirm ] = xlsWritePretty( sheets, names, fullfile)
%xlsWritePretty Write the excel sheets by going through csv
%   Sheets have sheets as seperate cells and are cells with the grid
%   names is a cell with all the names of the sheets


ind = regexp(fullfile,'[.]');
beg = fullfile(1:(ind(end)-1));

inddir = regexp(fullfile,'[/]');
begdir = fullfile(1:(inddir(end)));

if (~isunix())
    beg = strrep(beg, '/','\');
    fullfile = strrep(fullfile, '/','\');
    csvFile = [beg '.csv'];
end

realWd = strrep(pwd,'\','/');

for sh = 1: length(sheets)

    dms = size(sheets{sh});
    
    %Go through the different sheets if unix 
    %win can just use the same sheet, but now the names are the final sheet names
    
    if(isunix())
      csvList{sh} = strcat(begdir,names{sh});
      csvFile = csvList{sh};
    end


    %make the csv
    convertCell = cellfun(@(x)~(ischar(x) || isempty(x)),sheets{sh});
    sheets{sh}(convertCell) = cellfun(@(x)num2str(x, '%2.4f'),...
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



    %call VBscript or gnu ssconvert

    if (isunix())
      comd = ['ssconvert --import-encoding=stf_csvtab --export-type=Gnumeric_Excel:xlsx2 "' ...
       csvFile '" "' csvFile '"'];
    else
      comd = [realWd '/csv2xls.vbs "' fullfile '" ' int2str(sh) ' ' names{sh}];
    end

    [ranTest(sh)] = system(comd);
    if(logical(ranTest(sh)))
        error(['Error: command ' comd ' failed'])
    end

    csvLines = [];


end

%clean up
if (isunix())
  %combine all the produced csv files
  csvArgs = sprintf(' "%s" ',csvList{:});
  comd = ['ssconvert --merge-to "' fullfile '" ' csvArgs];
  [ranTest(sh+1)] = system(comd);
  for sh = 1: length(sheets)
    delete(csvList{sh})
  end
end
confirm = max(ranTest);

end
