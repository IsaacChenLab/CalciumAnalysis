function [BitMat,BitTimes]=MC_loadEvent(file)

global CAT;
global DATE;
global DATA_DIR;

BitMat=[]; BitTimes=[];

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  

fileName=MC_digFileName(file);

warning off all;
    
s=sprintf('load %s BitTimes BitMat;',fileName);
eval(s);

warning on;

return;