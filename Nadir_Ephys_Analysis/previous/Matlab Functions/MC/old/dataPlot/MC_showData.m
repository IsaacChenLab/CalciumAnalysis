function MC_showData(file,channels,Time,segment,notch60)
% function MC_showData(fileName,channels,Time,segment)
% shows channels (raw and highpass filtered) by segments of time, advance
% by mouse click.
% Time in seconds, so is segment

if nargin<4
    segment=1;
end

Fs=25e3;
if nargin<5
    notch60=0;
end

if notch60,
    [dataStream,totalSeconds,Hz]=MC_openFile(file);
    [B,A] = butter(2,[55 65]/(Hz/2),'stop');
end

[a1,t1]=MC_getElectrodes(file,channels,Time);
%o1=MC_filtElec(a1,Fs);
if notch60
    a1=filtfilt(B,A,a1);
end


sbn=[1 1; 2 1; 3 1; 2 2; 3 2; 3 2; 4 2; 4 2];
%seg=round(range(Time)/segment);
figure;
for t2=1:segment*Fs:length(t1)
    disptime=t2:t2+segment*Fs-1;
    for i=1:length(channels),
        subplot(sbn(length(channels),1),sbn(length(channels),2),i);
        hold on;
        plot(t1(disptime),a1(disptime,i),'k');
        %plot(t1(disptime),o1(disptime,i),'g');
        title(num2str(channels(i)));
        vv=[t1(disptime(1)) t1(disptime(end))];
        set(gca,'XLim',vv);
        %set(gca,'ylim',[min(a1(disptime,i))-range(a1(disptime,i))/5 max(a1(disptime,i))+range(a1(disptime,i))/5]);
        %set(gca,'XTick',[vv(1):0.1:vv(2)-eps vv(2)]);
    end
    %t=waitforbuttonpress;
end


return;

        

