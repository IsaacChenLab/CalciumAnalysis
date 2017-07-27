function [rawData,xtime,Hz,totalMs]=MC_loadElectrodes(file,ch,Time)

global CAT;
global DATE;
global DATA_DIR;

if nargin<3
    Time=nan;
end

global RAW_DATA_DIR;
if isempty(DATA_DIR)
    DATA_DIR='c:\data';
end
cd(DATA_DIR);
% if isempty(RAW_DATA_DIR)
%     RAW_DATA_DIR=DATA_DIR;
% end
% cd(RAW_DATA_DIR);
cd('PreparedFiles');

rawData=[];
xtime=[];
totalMs=0;
Hz=25e3;

sfl=sprintf('%s_%s_%d.mat',CAT,DATE,file);
if exist(sfl,'file');
    s=sprintf('load %s rawData_%d xtime totalMs Hz;',sfl,ch);
    eval(s);
    s=sprintf('rawData=rawData_%d; clear(''rawData_%d'');',ch,ch);
    eval(s);
    if ~isnan(Time)
        if xtime(end)>Time(end)
            inxs=find(xtime>=Time(1) & xtime<=Time(end));
            xtime=xtime(inxs);
            rawData=rawData(inxs);
        end
    end
else
    fileName=MC_fileName(file);
    [dataStream,totalMs,Hz]=MC_openFile(fileName);
    [rawData,xtime]=MC_getElectrodes(fileName,ch,Time);
end

return;

