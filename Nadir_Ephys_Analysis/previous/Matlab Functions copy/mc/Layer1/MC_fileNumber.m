function fileNumber=MC_fileNumber(fileName)
% function fileNumber=MC_fileNumber(fileName)
% gets a file name and retrieves it serial number

if isstr(fileName)
    fileNumber=MC_getNumber(fileName);
elseif iscell(fileName)
    for i=1:length(fileName)
        fileNumber(i)=MC_getNumber(fileName{i});
    end
elseif isnumeric(fileName)
    fileNumber=fileName;
end

return;

%%%%%%%%%%%%
function fileNumber=MC_getNumber(fileName)
global CAT;
global DATE;

if isempty(DATE),
    mark=(findstr(fileName,CAT)+length(CAT));
else
    mark=(findstr(fileName,CAT)+length(CAT)+length(DATE));    
end
mark=mark(end);
q=findstr(fileName,'.')-mark;
if q==0 %
    fileNumber=0;
else
    fileNumber=str2num(fileName(mark:mark+2));
end

return;
