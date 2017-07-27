function rez=BHV_Analyze (nm, num, method, param)

load (nm,'EVNTS','LCKS');
if ~exist ('EVNTS') || ~exist ('LCKS')
    error ('file does not contain required data');
end
switch (method)
    case 'hist'
        ref=EVNTS{param{1}};
        if param{2}==0
            evs=LCKS;
        else
            evs=EVNTS{param{2}};
        end
        
    case '%'
        ref=EVNTS{param{1}};
        rez=[];
        for i=1:length(ref)-num
            tm=ref(i:i+num);
            n=0;
            for j=1:length(tm)
                if ~isempty(find(EVNTS{1}>tm(j) & EVNTS{1}<tm(j)+5))
                    n=n+1;
                end
            end
            rez(end+1)=n/length(tm);
        end
    case 'disc'
        ev1=EVNTS{param{1}};
        ev2=EVNTS{param{2}};
        rez=[];
        e1=XCorr (ev1(end-num:end),LCKS,[-5 5],.1);
        p=find(ev2>ev1(end-num) & ev2<ev1(end));
        e2=XCorr (ev2(p),LCKS,[-5 5],.1);
        f=figure;
        hold on;
        plot(-5:0.1:5,e1);plot(-5:0.1:5,e2,'r');
        [x, ~, ~]=ginput(2);
        close (f);
        wnd=max(min(x),0.1);
        wnd(2)=max(x);
        wnd=round(wnd*10);
        h=waitbar(0,'Processing...');
        for i=1:length(ev1)-num
            waitbar(i/(length(ev1)-num));
            e1=XCorr (ev1(i:i+num),LCKS,[0 10],.1);
            p=find(ev2>ev1(i) & ev2<ev1(i+num));
            e2=XCorr (ev2(p),LCKS,[0 10],.1);
            e1=sum(e1(wnd(1):wnd(2)));
            e2=sum(e2(wnd(1):wnd(2)));
            rez(end+1)=(e1-e2)/(e1+e2);
        end
        close (h)
end
