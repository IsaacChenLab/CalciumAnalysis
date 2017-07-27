function [fileName,file]=MC_spShapeFileName(file)

global CAT;
global DATE;
global DATA_DIR;

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  
nfl=MC_multiFilesName(file);

fileName=sprintf('%s\\DigData\\SpikeShapes_%s.mat',DATA_DIR,nfl);

% if exist(fileName,'file')
%     fileName=[fileName ' -append'];
% end

return;
