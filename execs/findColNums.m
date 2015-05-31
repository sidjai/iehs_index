function [indVar indVal indUnit indField] = findColNums(raw)
%Find where the columns are for a input excel sheet
for h=1:size(raw,2)
    if (ischar(raw{1,h}) || ~isnan(raw{1,h}))
        if ~isempty(strfind(raw{1,h},'ble'))
            indVar=h;
        elseif ~isempty(strfind(raw{1,h},'alue'))
            indVal=h;
        elseif ~isempty(strfind(raw{1,h},'nit'))
            indUnit=h;
        elseif ~isempty(strfind(raw{1,h},'ield'))
            indField=h;
        end
    end
end
end
