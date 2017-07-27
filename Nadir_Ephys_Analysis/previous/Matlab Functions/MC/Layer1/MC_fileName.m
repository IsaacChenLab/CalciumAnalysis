function fileName=MC_fileName(file)
% function fileName=MC_fileName(file)
% builds a file name from file number

global DATA_DIR;
global CAT;
global DATE;
global RAW_DATA_DIR;
if isempty(DATA_DIR)
    DATA_DIR='c:\data';
end
if isempty(RAW_DATA_DIR)
    RAW_DATA_DIR=DATA_DIR;
end

if isnumeric(file)
    for i=1:length(file)
        if isempty(DATE)
            fileName{i}=sprintf('%s\\%s\\%s%03d.mcd',RAW_DATA_DIR,CAT,CAT,file(i));
        else
            fileName{i}=sprintf('%s\\%s\\%s\\%s%s%03d.mcd',RAW_DATA_DIR,CAT,DATE,CAT,DATE,file(i));
        end
    end
    if i==1
        fileName=fileName{i};
    end
else
    fileName=file;
end

% if ~exist(fileName,'file')
%     warning(['File ' fileName ' does not exist']);
% end

return;

    