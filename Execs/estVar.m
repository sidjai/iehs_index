function [flag] = estVar(chem, var)
%estVar say if you need to estimate this variable with some algorithm
%   0 has something in it, 1 = you should estimate
%   for the mixture (chem is a structure array): 1 means all the vars exist
%   and you should extimate the mixture property, 0 means something is missing


ln=length(chem);
flag=zeros(ln,1);
for i=1:ln
    if isfield(chem(i),var)
        val = chem(i).(var);
        if iscell(val)
            asdf= 0;
        else
            asdf=~isempty(find(isnan(val),1));
        end
        flag(i)= (isempty(val) || asdf);

    else
        flag(i)=1;
    end
end
flag=logical(flag);

end


