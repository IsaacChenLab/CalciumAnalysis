function [remMat,remTimes]=MC_loadREM(file)

global CAT;
global DATE;
global DATA_DIR;

remMat=[]; remTimes=[];

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  

fileName=MC_digFileName(file);

warning off all;
    
s=sprintf('load %s remTimes remMat;',fileName);
eval(s);

warning on;

return;