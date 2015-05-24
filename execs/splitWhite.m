function [val] = splitWhite(sin)
  sin = strtrim(sin);
  val = regexp(sin,'\s','split');
  numSet = cellfun(@(x)isempty(regexpi(x,'[a-df-z]','once')),val);
  val{numSet} = str2num(val{numSet});
  end
end