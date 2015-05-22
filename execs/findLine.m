function [val] = findLine(c , string)
val = find(~cellfun('isempty', strfind(c,string)));
end
