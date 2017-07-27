function MC_showData(fileName,channels,Time,segment)
% function MC_showData(fileName,channels,Time,segment)
% shows channels (raw and highpass filtered) by segments of time, advance
% by mouse click.
% Time in seconds, so is segment

Fs=25e3;

[a1,t1]=MC_getElectrodes(fileName,channels,Time);

sbn=[1 1; 2 1; 3 1; 2 2; 3 2; 3 2; 4 2; 4 2];

for t2=1:segment*Fs:length(t1)
    disptime=t2:t2+segment*Fs-1;
    for i=1:length(channels),
        subplot(sbn(length(channels),1),sbn(length(channels),2),i);
        hold on;
        plot(t1(disptime),a1(disptime,i),'b');
        vv=[t1(disptime(1)) t1(disptime(end))];
        set(gca,'XLim',vv);
        %set(gca,'XTick',[vv(1):0.1:vv(2)-eps vv(2)]);
    end
    t=waitforbuttonpress;
end


return;

        

