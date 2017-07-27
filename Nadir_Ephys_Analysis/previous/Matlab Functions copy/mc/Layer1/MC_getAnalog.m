function [anaData,xtime,Hz]=MC_getAnalog(fileName,anaNumbers,Time)
% function [anaData,xtime,Hz]=MC_getAnalog(fileName,anaNumbers,Time)
% inputs:
% fileName is an *.mcd file from MC_Rack
% anaNumbers: the original analog channel numbers in MC_Rack
% Time in seconds, otherwise gets the complete file
% outputs:
% anaData: a matrix where each column is an electrode
% xtime: a vector of the corresponding times in seconds
% Hz: sampling rate

if nargin<3
    Time=nan;
end

[anaData, xtime, Hz]=MC_getAnalogChannels(fileName,'Analog Raw Data',anaNumbers,Time);

return;

