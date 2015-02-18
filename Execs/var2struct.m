function s = var2struct(varargin)
  names = arrayfun(@inputname,1:nargin,'UniformOutput',false);
  s = cell2struct(varargin,names,2);
end