function [ out ] = xlsReadPretty(fullfile)
%xlsReadPretty Read xls documents by converting to csv
%   Detailed explanation goes here

realWd = strrep(pwd,'\','/');
comd = [realWd '/xls2csv "' fullfile '"'];
[ranTest completeTest] = system(comd);
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
        commas = regexp(parse{row},',');
        col = 2;
        bef = commas(1);
        if bef~=1
            out{sh}{row,1} = parse{row}(1:bef-1);
        end
        
        for ca = commas(2:end)
            %if the commas are next to each other just add col
            
            if (bef+1 ~= ca)
               out{sh}{row,col} = parse{row}(bef+1:ca-1);
               
            end
            col=col+1;
            bef = ca;
        end
        
    end
    
    
    delete(csvFile)

    sh = sh+1;
    csvFile = [beg num2str(sh) '.csv'];
    
end

