function [rms,threshold,outliers,whatevents,evinterv]=MC_cutSpikesGetParam(rawData,Times,Fs,digData)
% function [waveForms,spikeTimes]=MC_cutSpikes(rawData,Times,Fs,Threshold)
% filters the rawData with a high pass filter, finds thresholds crossings,
% and extracts -1 till +2 ms as the spikes waveforms.
% inputs:
% rawData is a vector of voltages
% Times is the corresponding time in second
% Fs assumed to be 25e3, unless supplied
% procedure:
% 1. filters high pass
% 2. shows data and asks wether RMS or not should be used
% 3. asks for outliers threshold i.e. a y-axis value that if corssed, the
% 'fake spike' will not counted.
% 4. asks for threshold, i.e. a y-axis value that if corossed the event
% will be counted as spike.
% 5. extract the spikes wave forms, -2ms to +3ms.
% 6. realign them on minimum voltage and takes -1ms to +2ms

if nargin<3
    digData=[];
end
if nargin<4,
    Fs=25e3;
end
THRESHOLD=6;
rms='n';
MIN_STD=6;
OUTLIERS=15;
c={'b','r','g','k','c','m','y','b'};
whatevents=[];
evinterv=[];

fltData=MC_spHighPass(rawData,Fs);
clear rawData;

if ~isempty(digData) & ~isnan(digData)
    evs=unique(digData(:)); evs(evs==0)=[];
else
    evs=[];
end

[rmsData,rmsTime,rmsSamples]=MC_spRMS(fltData,Times,Fs);
scrsz = get(0,'ScreenSize');
f1=figure('Position',[round(scrsz(3)/10) round(scrsz(4)/10) round(8*scrsz(3)/10) round(8*scrsz(4)/10)]);
b1=0; b2=0;
for i=1:size(rmsData,2)
    subplot(1,2,2); hold on;
    plot(rmsTime,rmsData(:,i)+b1);
    b1=b1+30*std(rmsData(:,i));
    set(gca,'xlim',[rmsTime(1) rmsTime(end)]);
    set(gca,'ylim',[0 b1]);
    subplot(1,2,1); hold on;
    plot(Times,fltData(:,i)+b2);
    b2=b2+30*std(fltData(:,i));
    set(gca,'xlim',[Times(1) Times(end)]);
    set(gca,'ylim',[ -b2/(i+1) b2]);
    % plot events
    for j=1:length(evs)
        x=find(digData(:,i)==evs(j));
        plot(Times(x),b2*ones(length(x),1),'v','color',c{j});
    end
end

rms=questdlg('','','RMS','Amplitude','Cancel','Amplitude');
if strcmp(rms,'Cancel')
    close all;
    error('no selection');
elseif strcmp(rms,'RMS')
    rms=1;
else
    rms=0;
end

num_files=size(fltData,2); for i=1:num_files, str{i}=num2str(i); end;
[fl,ok] = listdlg('PromptString','Based on file:','SelectionMode','single','ListString',str);
if ~ok
    error('stop');
else
    file_n=str2num(str{fl});
end

close(f1);
if rms,
    cData=rmsData(:,file_n);
    clear rmsData;
    xaxis=rmsTime;
else
    cData=fltData(:,file_n);
    clear fltData;
    xaxis=Times;
end
dData=digData(:,file_n);

c_evs=[];
ccData=cData; cdData=dData; cxaxis=xaxis;
f1=figure('Position',[round(scrsz(3)/10) round(scrsz(4)/10) round(8*scrsz(3)/10) round(8*scrsz(4)/10)]);
while 1
    cla; hold on;
    plot(xaxis,cData);
    set(gca,'xlim',[xaxis(1) xaxis(end)]);
    if rms
        set(gca,'ylim',[0 max(mean(cData)+30*std(cData),max(cData))]);
    else
        set(gca,'ylim',[min(min(cData),mean(cData)-40*std(cData)) max(max(cData),mean(cData)+40*std(cData))]);   
    end
    if isempty(evs)
        break;
    end
    str=[];
    ymin=get(gca,'ylim'); ymin=ymin(1);
    for j=1:length(evs)
        x=find(dData==evs(j));
        plot(xaxis(x),ymin*ones(length(x),1),'^','color',c{j},'markersize',9);
        str{j}=[c{j} ' - ' num2str(evs(j))];
    end
    str{end+1}='Return';
    if MC_accept 
        break;
    end
    [c_evs,ok] = listdlg('PromptString','Events?','SelectionMode','mutiple','ListString',str);
    if ~ok
        error('stop');
    end
    if ~strcmp(str(c_evs),'Return')
        while 1
            s=inputdlg({'from','to'},'Interval to dismiss (ms)?',ones(2,1),{'1','100'});
            interv1=str2num(s{1});
            interv2=str2num(s{2});
            if isnumeric(interv1) & isnumeric(interv2)
                break;
            end
        end
        for i=1:length(c_evs)
            x=find(whatevents==evs(c_evs(i)));
            if isempty(x)
                whatevents(end+1)=evs(c_evs(i));
                evinterv(end+1,:)=[interv1 interv2];
            else
                evinterv(x,:)=[interv1 interv2];
            end
        end
        evinterv=(Fs/1e3)*evinterv;
        [cData]=MC_spRemoveAroundEvents(cData,dData,xaxis,whatevents,evinterv);
    else
        whatevents=[]; evinterv=[];
        cData=ccData; dData=cdData; xaxis=cxaxis;
    end
       
end
if rms
    [x,y]=ginput(1);   
else
    [x,y]=ginput(2);
end
if ~isempty(y),
    outliers=sort(y);
else
    if rms
        outliers=mean(cData)+OUTLIERS*std(cData);
    else
        outliers=[mean(cData)-OUTLIERS*std(cData) mean(cData)+OUTLIERS*std(cData)];
    end
end
ll=get(gca,'ylim');
if rms, 
    set(gca,'ylim',[0 y]); 
else
    set(gca,'ylim',outliers); 
end
[x,y]=ginput(1);
if ~isempty(y)
    threshold=y;
else
    if rms
        threshold=mean(cData)+THRESHOLD*std(cData);
    else
        threshold=mean(cData)-THRESHOLD*std(cData);
    end
end
close(f1);


return;
