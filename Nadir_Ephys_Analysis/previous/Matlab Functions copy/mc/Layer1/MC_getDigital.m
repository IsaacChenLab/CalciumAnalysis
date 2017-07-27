function [digData,xtime,Hz]=MC_getDigital(fileName,Time)
% function digTimes=MC_getDigitalTimes(fileName,Time)
% inputs:
% fileName is an *.mcd file from MC_Rack
% Time in mili-seconds, otherwise gets the complete file
% outputs:
% digTimes: a cell array where each cell is a vector with time in mili-seconds of onset of
% the bit.

[dataStream,totalMS,Hz]=MC_openFile(fileName);

if nargin<2 | isnan(Time) | strcmp(Time,'all')
    Time=[0 totalMS];
end

ce = nextdata(dataStream,'streamname','Digital Data','startend',Time);

digData=ce.data(1,:)';

ll=size(digData,1);
xtime=Time(1)+(1e3/Hz)*[1:ll]';

return;
