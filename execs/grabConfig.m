function [ input ] = grabConfig( config )
%grabConfig Sets the configuration for the run using the text file.
%   If a file name has a : then it will be taken as a full location,
%   otherwise it will be added to the direc
%   New vars can be added by making a new line with 'var = val' formulation
%   without quotes or semicolons.
%   Also grabs the global variables and templates that are useful

fid = fopen(config);
C= textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

findLine = @(c , string) (find(~cellfun('isempty', strfind(c,string))));
raw = C{1};
useSet = findLine(raw,'=')';
comSet = findLine(raw,'%');

useSet = useSet(~ismember(useSet,comSet));

dirSet = findLine(raw,'direc');
direc = raw{dirSet}(strfind(raw{dirSet},'=')+1:end);
input.direc = direc;
useSet = useSet(~ismember(useSet, dirSet));

for n = useSet
    [ref, val] = strtok(raw{n}, '=');
    ref = strtrim(ref);
    val = strtrim(val(2:end));
    
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
input.AppendixName = [input.IEHSdir '\Templates\IO_notes.xlsx'];
input.MasterChemListName = [input.IEHSdir '\GlobalVars\masterChemList.txt'];
input.ClassDesName = [input.IEHSdir '\GlobalVars\classdesignations.xlsx'];

% load the appendix and templates
for sh=2:4
    [txt,txt,junk]=xlsread(input.AppendixName,sh,'A3:H30');
    IOnotes{sh}=txt;
    for pa=1:size(txt,1)
        nte = find(~cellfun('isempty',txt(pa,:)),1,'last');
        if isempty(nte)
            lasts (pa,sh) = 0;
        else
            lasts(pa,sh) = nte;
        end
    end    
end
IOnotes{1}=lasts;
input.IOnotes = IOnotes;

end

