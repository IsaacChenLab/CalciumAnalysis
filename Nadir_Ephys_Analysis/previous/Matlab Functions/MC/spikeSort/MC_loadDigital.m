function [digData,xtime,Hz,totalMs]=MC_loadDigital(file,Time)
% function digTimes=MC_getDigitalTimes(fileName,Time)
% inputs:
% fileName is an *.mcd file from MC_Rack
% Time in mili-seconds, otherwise gets the complete file
% outputs:
% digTimes: a cell array where each cell is a vector with time in mili-seconds of onset of
% the bit.

global CAT;
global DATE;
global DATA_DIR;

if nargin<2
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

digData=[];
xtime=[];
totalMs=0;
Hz=25e3;

sfl=sprintf('%s_%s_%d.mat',CAT,DATE,file);
if exist(sfl,'file');
    s=sprintf('load %s digData xtime totalMs Hz;',sfl);
    if exist('digData','var')
        eval(s);
        if ~isnan(Time)
            if xtime(end)>Time(end)
                inxs=find(xtime>=Time(1) & xtime<=Time(end));
                xtime=xtime(inxs);
                digData=digData(inxs);
            end
        end
    else
        fileName=MC_fileName(file);
        [dataStream,totalMs,Hz]=MC_openFile(fileName);
        [digData,xtime]=MC_getDigital(fileName,Time);
    end
else
    fileName=MC_fileName(file);
    [dataStream,totalMs,Hz]=MC_openFile(fileName);
    [digData,xtime]=MC_getDigital(fileName,Time);
end

return;
