function cdc

global DATA_DIR;
global CAT;
global DATE;

cd(DATA_DIR);
cd(CAT);
if ~isempty(DATE)
    cd(DATE);
end

return;

