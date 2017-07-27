function [digTimes,digInx]=MC_getDigitalTimes(fileName,bits,Time)
% function digTimes=MC_getDigitalTimes(fileName,bits,Time)
% inputs:
% fileName is an *.mcd file from MC_Rack
% bits: the digital bits numbers in MC_Rack
% Time in mili-seconds, otherwise gets the complete file
% outputs:
% digTimes: a cell array where each cell is a vector with time in mili-seconds of onset of
% the bit.

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

if nargin<2
    bits=1;
end

digData=[];
for i=1:periods,
    ce = nextdata(dataStream,'streamname','Digital Data','startend',Start+[0 CHUNKS]+(i-1)*CHUNKS);
     ddata=ce.data(bits,:)';
    digData=[digData; ddata];
end
ce = nextdata(dataStream,'streamname','Digital Data','startend',Start+[periods*CHUNKS periods*CHUNKS+last]);
ddata=ce.data(bits,:)';
digData=[digData; ddata];

for i=1:length(bits),
    x=find(digData(:,i)>0);
    x=runs(x);
    digTimes{i}=(1e3*x(1,:)/(Hz))';
    digInx{i}=x(1,:)';
end

if length(bits)==1
    digTimes=digTimes{1};
    digInx=digInx{1};
end

return;
