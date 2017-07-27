function nch=MC_NumElecInFile(fileName)
% function nch=MC_NumElecInFile(fileName)
% retrieves the number of recorded electrodes in the file

[dataStream,totalSeconds,Sampling_Hz]=MC_openFile(fileName);
datastrms=getfield(dataStream,'StreamNames');
datainfos=getfield(dataStream,'StreamInfo');
nch=getfield(dataStream,'NChannels2');
streamName='Electrode Raw Data';
i=strmatch(streamName,datastrms,'exact');
nch=nch(i);

return;

