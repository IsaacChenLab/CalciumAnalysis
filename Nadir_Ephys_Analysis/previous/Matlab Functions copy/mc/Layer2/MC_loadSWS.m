function [swsMat,swsTimes]=MC_loadSWS(file)

global CAT;
global DATE;
global DATA_DIR;

swsMat=[]; swsTimes=[];

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  

fileName=MC_digFileName(file);

warning off all;

s=sprintf('load %s swsTimes swsMat;',fileName);
eval(s);

warning on;


return;