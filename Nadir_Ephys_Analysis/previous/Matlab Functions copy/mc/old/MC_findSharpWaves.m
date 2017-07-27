function [Times,Mat]=MC_findSharpWaves(file,ch,Time)

if nargin<3
    Time=nan;
end

C=5000; % chunks of 5
Total_ms=0;
Times=[];
for i=1:length(file),
    
    file=MC_fileNumber(file);
    [dataStream,totalMs,Hz]=MC_openFile(file(i));

    [rawData,xtime]=MC_getElectrodesInMs(file(i),ch,Time);
    rawData=MC_notch60(rawData,Hz);

%     [B,A] = butter(2,[300]/(Hz/2));
%     rawData=filtfilt(B,A,rawData);

    periods=ceil(xtime(end)/C);
    last=mod(xtime(end),C);

    d=0; st=[];
    f1=figure;
    for i=1:periods,
        if last>0 & i==periods
            tm=[0+[0 last]+(i-1)*C];
        else
            tm=[0+[0 C]+(i-1)*C];
        end
        curSegInx=find(xtime>=tm(1) & xtime<tm(2));
        curtime=xtime(curSegInx);
        curdata=rawData(curSegInx);
        plot(curtime,curdata,'k'); 
        hold on;
        while (1)
            [x,y]=ginput(1);
            if isempty(x)
                break;
            end
            [cur_wave,inx,inxs]=MC_findSharpNeigh(curtime,curdata,x);
            plot(curtime(inxs),curdata(inxs),'r');
            plot(curtime(inx),curdata(inx),'ro');
            s=questdlg('Accept')
            if strcmp(s,'Yes')
                d=d+1;
                st(d)=cur_wave;
            else
                plot(curtime(inx),curdata(inx),'o','color',[1 1 1]);
                plot(curtime(inxs),curdata(inxs),'k');
            end    
        end
        hold off;
    end
    close(f1);
    
    Times=[Times Total_ms+st];
    Total_ms=Total_ms+totalMs;

    
end


Mat=sparse(Total_ms,1);
Mat(Times,1)=1;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Time,inx,inxs]=MC_findSharpNeigh(tt,et,x)

R=150; % ms

% search for minimum within R ms

minus_R=x-R;
if minus_R<tt(1),
    minus_R=tt(1);
end
plus_R=x+R;
if plus_R>tt(end),
    plus_R=tt(end);
end
inxs=find(tt>=minus_R & tt<plus_R);
[mn,inx]=min(et(inxs));
inx=inx+inxs(1)-1;

Time=tt(inx);

return;





