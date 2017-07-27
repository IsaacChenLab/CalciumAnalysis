function rez=DA_Calibration(fl, con, method)

load (fl)
switch lower(method)
    
    case 'val'
        figure;
        plot (SWEEP_DATA(1:1200,:),'k');
        title ('Choose the position')
        [val, ~]=ginput(1);
        close;
        data=SWEEP_DATA(round(val),:);
        clear SWEEP_DATA;
end

figure;

rep=1;
col={'r'; 'g'; 'b'; 'm'; 'y'};
while rep
    hold off
    plot (data,'k')
    hold on
    ran=[];ran(length(con),2)=0;
    for i=1:length(con)
        title (['Choose period for [DA] = ' num2str(con(i))]);
        [xs, ~]=ginput(2);
        ran(i,1)=round(min(xs));
        ran(i,2)=round(max(xs));
        plot (ran(i,1):ran(i,2),data(ran(i,1):ran(i,2)),col{i});
    end
    a=questdlg ('Are the periods acceptable?','Question','Yes','No','No');
    if strcmp(a,'Yes')
        rep=0;
    end
end
close
rez(1:length(con),1)=con;
for i=1:length(con)
    rez(i,2)=mean(data(ran(i,1):ran(i,2)))
end
figure;plot(rez(:,1),rez(:,2),'kx');