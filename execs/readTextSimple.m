function [cellstr] =  readTextSimple(path, head)


  if(~exist(path,'file'))
    error(['File: ' path ' does not exist']);
  end

  fid = fopen(path);
  cellstr = textscan(fid, '%s', 'Delimiter', '\n','HeaderLines',head);
  cellstr = strtrim(cellstr{1});

  fclose(fid);
end
