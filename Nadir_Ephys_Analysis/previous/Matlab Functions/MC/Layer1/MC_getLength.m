function ms=MC_getLength(file)

ms=0; 
for i=1:length(file)
    % get raw data
    [rawData,xtime,Hz,totalMs]=MC_loadElectrodes(file(i),-1);
    ms=ms+totalMs;
end

return;