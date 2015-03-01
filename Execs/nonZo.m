function [ val ] = nonZo( vin )
%nonZo Return the nonzero values
%   if they are all zero, output 0
ind = logical(vin);
val = vin(ind);
if(isempty(val))
    val = 0;
end
end

