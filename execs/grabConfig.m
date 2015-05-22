function [ input ] = grabConfig( config )
%grabConfig Sets the configuration for the run using the text file.
%   If a file name has a : then it will be taken as a full location,
%   otherwise it will be added to the direc
%   New vars can be added by making a new line with 'var = val' formulation
%   without quotes or semicolons.
%   Also grabs the global variables and templates that are useful

fid = fopen(config);
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

raw = C{1};
useSet = findLine(raw,'=')';
comSet = findLine(raw,'%');

useSet = useSet(~ismember(useSet,comSet));
dirSet = findLine(raw,'direc');
direc = regexprep(raw{dirSet}, '\\', '/')
direc = direc(strfind(direc,'=')+1:end);
input.direc = direc;
useSet = useSet(~ismember(useSet, dirSet));

for n = useSet
    [ref, val] = strtok(raw{n}, '=');
    ref = strtrim(ref);
    val = strtrim(val(2:end));
    val = regexprep(val, '\\', '/');
    
    if (isempty(regexpi(val,'[a-z]','once')))
        val = str2double(val);
    elseif (isempty(strfind(val,':')))
        val = [direc val];
    end
    
    input.(ref) = val;
    
    
    
end

%Gets location for the Global variable files
workdir = pwd;
input.IEHSdir = fileparts(workdir);
input.AppendixName = [input.IEHSdir '/templates/IO_notes.xlsx'];
input.MasterChemListName = [input.IEHSdir '/global_vars/masterChemList.txt'];
input.ClassDesName = [input.IEHSdir '/global_vars/classdesignations.xlsx'];

% load the appendix and templates
app = xlsReadPretty(input.AppendixName, 999, 'header');
for sh=2:4
    for pa=1:size(app,1)
        nte = find(~cellfun('isempty',app(pa,:)),1,'last');
        if isempty(nte)
            lasts (pa,sh) = 0;
        else
            lasts(pa,sh) = nte;
        end
    end    
end
input.IOnotes = {lasts, app{2:4}};

end

