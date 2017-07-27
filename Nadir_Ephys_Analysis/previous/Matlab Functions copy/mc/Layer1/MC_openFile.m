function [dataStream,totalMs,Sampling_Hz]=MC_openFile(file)
% function [dataStream,totalSeconds,Sampling_Hz]=MC_openFile(file)
% open an MC_Rack file, if fileName is a number opens the file
% c[CAT#][file#].mcd
% returns the data stream, total recorded time in seconds, and the sampling
% frequency in Hertz


fileName=MC_fileName(file);

% open file with MCrack interface software
a=datastrm(fileName);

totalMs=round(1e3*24*3600*(getfield(a,'recordingStopDate')-getfield(a,'recordingdate')));
totalMs=totalMs-1e3;

mst=getfield(a,'MicrosecondsPerTick');
Sampling_Hz= 1e6 / mst;

dataStream=a;

return;