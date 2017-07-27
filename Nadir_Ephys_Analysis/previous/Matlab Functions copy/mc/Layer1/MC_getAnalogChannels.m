function [analogData, xtime, Hz]=MC_getAnalogChannels(fileName,streamName,channelNumbers,Time)
% function analogData=MC_getAnalogChannels(fileName,streamName,channelNumbers,Time)
% inputs:
% fileName is an *.mcd file from MC_Rack
% exampls for streamName: 'Electrode raw data',' Filtered Data 1'
% channelNumbers: the original (electrode) numbers in MC_Rack
% Time in seconds, otherwise gets the complete file
% outputs:
% a matrix where each column is a channel

CHUNKS=10000; % chunks of 10000 ms

[dataStream,totalMS,Hz]=MC_openFile(fileName);

if nargin<4 | isnan(Time) | strcmp(Time,'all')
    Start=0;
    TimeS=[0 totalMS];
else
    TimeS=Time;
    Start=Time(1);
    totalMS=Time(2)-Time(1);
end

periods=floor(totalMS/CHUNKS);
last=mod(totalMS,CHUNKS);

analogData=[];
for i=1:periods,
    ce = nextdata(dataStream,'streamname',streamName,'originorder','on','startend',Start+[0 CHUNKS]+(i-1)*CHUNKS);
    [vdata, tvals] = ad2muvolt(dataStream, ce.data(channelNumbers,:), streamName);
    vdata=vdata';
    analogData=[analogData; vdata];
end
ce = nextdata(dataStream,'streamname',streamName,'originorder','on','startend',Start+[periods*CHUNKS periods*CHUNKS+last]);
[vdata, tvals] = ad2muvolt(dataStream, ce.data(channelNumbers,:), streamName);
vdata=vdata';
analogData=[analogData; vdata];
clear ce tvals vdata;
ll=size(analogData,1);
xtime=Start+(1e3/Hz)*[1:ll]';

return;
