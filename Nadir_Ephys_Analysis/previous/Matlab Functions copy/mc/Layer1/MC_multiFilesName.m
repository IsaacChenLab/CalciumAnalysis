function nfl=MC_MultiFilseName(file)

global DATA_DIR;
global CAT;
global DATE;

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  
if length(file)==1 | ~isempty(find(diff(file)>1))
    nfl=sprintf('%d',file);
else
    nfl=sprintf('%d-%d',file(1),file(end));
end

nfl=sprintf('%s_%s_%s',CAT,DATE,nfl);

return;

    