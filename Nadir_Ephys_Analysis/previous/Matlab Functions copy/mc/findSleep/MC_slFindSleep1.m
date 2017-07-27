function [swsMat,swsTimes,remMat,remTimes,FileTimes]=MC_slFindSleep(file,ch,eog,RMS_MS)

if nargin<3
    eog=nan; eogData=nan;
end

if nargin<4
    RMS_MS=2000; % ms for RMS
end

SEC_CHUNKS = 60e3; %ms for presentation in chunks

%try,
[swsMat,swsTimes,remMat,remTimes,FileTimes]=deal(nan);
[swsTimes,remTimes]=deal([]);
f1=figure('Units','normalized','position',[0.15 0.1 0.7 0.8]);
Time=[]; rawData=[]; eogData=[]; events=[];
Total_ms=0;
% over files
for i=1:length(file)
    % get raw data
    [rd,totalMs,Hz,T]=MC_loadLFP(file(i),ch);
    rawData=[rawData; rd];
    Time=[Time; T+Total_ms];
    if ~isnan(eog)
        [ed]=MC_loadLFP(file(i),eog);
    else
        ed=zeros(size(rd));
    end
    clear rd T;
    eogData=[eogData; ed];
    clear ed;
    [ev,tt,dHz]=MC_loadDigital(file(i)); clear tt;
    ev=MC_bin(ev,dHz/Hz);
    events=[events; ev];

    Total_ms=Total_ms+totalMs;
    FileTimes(i)=Total_ms;
end

rmsData=MC_slRMS(rawData,RMS_MS*Hz/1e3);
eogData=MC_BandPass(eogData,Hz,[1 100]);

rmslimits=[min(rmsData) max(rmsData(100:end-100)) mean(rmsData) std(rmsData)];
rawlimits=[min(rawData) max(rawData(100:end-100)) mean(rawData) std(rawData)];
eoglimits=[min(eogData) max(eogData(100:end-100)) mean(eogData) std(eogData)];

%[events]=MC_loadEvents(file);
%events=(sum(events,2)>0);

CHUNKS=SEC_CHUNKS * Hz / 1e3;
periods=floor(size(rawData,1)/CHUNKS);
last=mod(size(rawData,1),CHUNKS);
for i=1:periods,
    inx=[1:CHUNKS]+(i-1)*CHUNKS;
    [st,rt]=MC_slFindTimes(Time(inx),rmsData(inx),rmslimits,rawData(inx),...
        rawlimits,eogData(inx),eoglimits,events(inx),f1);   
    swsTimes=[swsTimes; st];
    remTimes=[remTimes; rt];
end  
inxf=[periods*CHUNKS+1 : periods*CHUNKS+last];
filler=nan(length(inx)-length(inxf),1); 
[st,rt]=MC_slFindTimes(Time(inx)+SEC_CHUNKS,[rmsData(inxf); filler],rmslimits,...
    [rawData(inxf); filler],rawlimits,[eogData(inxf); filler],eoglimits,[events(inxf); filler],f1);    
swsTimes=[swsTimes; st];
remTimes=[remTimes; rt];

if ~MC_accept('save?')
    [swsMat,swsTimes,remMat,remTimes,FileTimes]=deal(nan);
    close(f1);
    return;
end

clear rmsData rawData eogData Time;

close(f1);

% merge close times for edges of presentation, based on minimal Hz
if ~isempty(swsTimes)
    df=swsTimes(2:end,1)-swsTimes(1:end-1,2);
    x=find(df<=1.5*1e3/Hz);
    swsTimes(x,2)=swsTimes(x+1,2);
    swsTimes(x+1,:)=[];
end
if ~isempty(remTimes)
    df=remTimes(2:end,1)-remTimes(1:end-1,2);
    x=find(df<=1.5*1e3/Hz);
    remTimes(x,2)=remTimes(x+1,2);
    remTimes(x+1,:)=[];
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


% catch,
%     warning(lasterr);
%     [swsMat,swsTimes,remMat,remTimes,FileTimes]=deal(nan);
%     return;
% end
% 
return;


