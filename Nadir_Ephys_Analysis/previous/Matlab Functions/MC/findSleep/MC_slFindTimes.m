function [st,rt]=MC_slFindTimes(time,rmsData,rmslimits,rawData,rawlimits,eogData,eoglimits,events,f1);

if nargin<8
    f1=figure;
end

st=[]; rt=[];

if ~isempty(rawData)
    subplot(3,1,1); cla; hold on;
    plot(time,rawData,'b');
    set(gca,'xlim',[time(1) time(end)]);
    set(gca,'ylim',[(rawlimits(3)-7*rawlimits(4)) (rawlimits(3)+7*rawlimits(4))]);
    plot(time,rawlimits(3)*ones(size(rawData,1),1),'k');
    plot(time,(rawlimits(3)+2*rawlimits(4))*ones(size(rawData,1),1),'--k');
    plot(time,(rawlimits(3)-2*rawlimits(4))*ones(size(rawData,1),1),'--k');
    if ~isempty(events)
        a=find(events); a=MC_runs(a); a=a(1,:);
        plot(time(a),(rawlimits(3)-6.5*rawlimits(4))*ones(length(a),1),'r^','markersize',6);
    end
end

subplot(3,1,2); cla; hold on;
d1=rmsData; clear rmsData;
plot(time,d1,'b');
plot(time,rmslimits(3)*ones(size(d1,1),1),'k');
plot(time,(rmslimits(3)+2*rmslimits(4))*ones(size(d1,1),1),'--k');
set(gca,'xlim',[time(1) time(end)]);
set(gca,'ylim',[rmslimits(1) (rmslimits(3)+5*rmslimits(4))]);

shift=2e5; top=-1e5;

subplot(3,1,3); cla; hold on;
d2=eogData-shift; clear eogData;
eoglimits=eoglimits-[shift shift shift eps];
plot(time,d2,'b');
set(gca,'xlim',[time(1) time(end)]);
set(gca,'ylim',[(eoglimits(3)-20*eoglimits(4)) (eoglimits(3)+20*eoglimits(4))]);
plot(time,eoglimits(3)*ones(size(d2,1),1),'k');
plot(time,(eoglimits(3)+2*eoglimits(4))*ones(size(d2,1),1),'--k');
plot(time,(eoglimits(3)-2*eoglimits(4))*ones(size(d2,1),1),'--k');

ind=1;
while (ind)
    figure(f1);
    [x,y]=ginput(2);
    if isempty(x) | length(x)==1
        break;
    end
    if x(1)<time(1)
        x(1)=time(1);
    end
    if x(2)>time(end)
        x(2)=time(end);
    end
    inx=[find(time>=x(1),1) : find(time<=x(2),1,'last')];
    if y(1)<top, % eog
        subplot(3,1,3);
        plot(time(inx),d2(inx),'r');
        if MC_accept
            rt=[rt; x(1) x(2)];
        else
            plot(time(inx),d2(inx),'b');
        end        
    else
        subplot(3,1,2);
        plot(time(inx),d1(inx),'r');
        f2=figure('Units','normalized','position',[0.05 0.05 0.4 0.3]);
        plot(time(inx),rawData(inx)); 
        set(gca,'xlim',[time(inx(1)) time(inx(end))]);
        if MC_accept
            st=[st; x(1) x(2)];
        else
            figure(f1);
            subplot(3,1,2);
            plot(time(inx),d1(inx),'b');
        end
        close(f2);
    end
    
end

if nargin<7
    close(f1);
end

return;
