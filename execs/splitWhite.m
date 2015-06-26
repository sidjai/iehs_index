function [val] = splitWhite(sin)
  sin = strtrim(sin);
  val = regexp(sin,'\s','split');
  numSet = cellfun(@(x)isempty(regexpi(x,'[a-df-z]','once')),val);
  for si = find(numSet)
      val{si} = str2num(val{si});
  end
end