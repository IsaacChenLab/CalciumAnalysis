function [swsMat,swsTimes,remMat,remTimes,FileTimes]=MC_slFindSleep(file,ch,eog)

if nargin<3
    eog=nan; eogData=nan;
end


SEC_CHUNKS = 30e3; %ms

%try,

[swsTimes,remTimes]=deal([]);
f1=figure('Units','normalized','position',[0.15 0.2 0.7 0.6]);

Total_ms=0;
% over files
for i=1:length(file)
    % get raw data
    [rawData,totalMs,Hz,Time]=MC_loadLFP(file(i),ch);
    rmsData=MC_slRMS(rawData,1000*Hz/1e3);
    mm=max(rmsData); mn=min(rmsData); limits=[mn mm];
    mmn=mean(rmsData); tpp=2*std(rmsData);
    if ~isnan(eog)
        [eogData]=MC_loadLFP(file(i),eog);
    else
        eogData=zeros(size(rawData));
    end
    CHUNKS=SEC_CHUNKS * Hz / 1e3;
    periods=floor(size(rawData,1)/CHUNKS);
    last=mod(size(rawData,1),CHUNKS);
    for i=1:periods,
        inx=[1:CHUNKS]+(i-1)*CHUNKS;
        [st,rt]=MC_slFindTimes(rmsData(inx),rawData(inx),eogData(inx),Time(inx),limits,[mmn tpp],f1);    
        swsTimes=[swsTimes; st+Total_ms];
        remTimes=[remTimes; rt+Total_ms];
    end  
    inxf=[periods*CHUNKS+1 : periods*CHUNKS+last];
    filler=nan(length(inx)-length(inxf),1); 
    [st,rt]=MC_slFindTimes([rmsData(inxf); filler],[rawData(inxf); filler],[eogData(inxf); filler],Time(inx)+SEC_CHUNKS,limits,[mmn tpp],f1);    
    swsTimes=[swsTimes; st+Total_ms];
    remTimes=[remTimes; rt+Total_ms];
    
    Total_ms=Total_ms+totalMs;
    FileTimes(i)=Total_ms;
    clear rawData eogData;
end


% final sparse matrix of ones where the spike has occured
swsMat=sparse(Total_ms,1);
for i=1:size(swsTimes,1)
    swsMat(round(swsTimes(i,1)):round(swsTimes(i,2)))=1;
end
remMat=sparse(Total_ms,1);
for i=1:size(remTimes,1)
    remMat(round(remTimes(i,1)):round(remTimes(i,2)))=1;    
end

close(f1);

% catch,
%     warning(lasterr);
%     swsMat=nan;
%     swsTimes=nan;
%     remMat=nan;
%     remTimes=nan;
%     FileTimes=nan;
%     return;
% end
% 
return;


