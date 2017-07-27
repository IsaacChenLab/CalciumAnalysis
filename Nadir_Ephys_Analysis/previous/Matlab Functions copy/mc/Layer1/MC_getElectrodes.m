function [elecData,xtime,Hz]=MC_getElectrodes(fileName,elecNumbers,Time)
% function [elecData,xtime,Hz]=MC_getElectrodes(fileName,elecNumbers,Time)
% inputs:
% fileName is an *.mcd file from MC_Rack
% elecNumbers: the original electrode numbers in MC_Rack
% Time in seconds, otherwise gets the complete file
% outputs:
% elecData: a matrix where each column is an electrode
% xtime: a vector of the corresponding times in seconds

if nargin<3
    Time=nan;
end

[elecData, xtime, Hz]=MC_getAnalogChannels(fileName,'Electrode Raw Data',elecNumbers,Time);

return;

