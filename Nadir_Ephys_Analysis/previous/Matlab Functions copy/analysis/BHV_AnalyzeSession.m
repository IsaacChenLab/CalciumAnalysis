function rez=BHV_AnalyzeSession (fn)

try
    load (fn);
catch
    error ('file could not be loaded');
end
if ~(exist('EVENTS','var') && exist ('LICKS','var'))
    error ('the loaded session does not contain appropriate variables');
end
ev1=MC_smooth(XCorr (EVENTS{1},LICKS,[-5 10],.1)',0.5);
ev2=MC_smooth(XCorr (EVENTS{2},LICKS,[-5 10],.1)',0.5);
ev3=MC_smooth(XCorr (EVENTS{3},LICKS,[-5 10],.1)',0.5);
rt=Rate(LICKS,10)/10;
mx=max(max(max(ev1),max(ev2)),max(ev3(:)));

figure
subplot (2,2,1)
plot (-5:.1:10,ev1);
axis ([-5 10 0 mx+0.1]);
title ('EVENT 1');
subplot (2,2,2)
plot (-5:.1:10,ev2);
axis ([-5 10 0 mx+0.1]);
title ('EVENT 2');
subplot (2,2,3)
plot (-5:.1:10,ev3);
axis ([-5 10 0 mx+0.1]);
title ('EVENT 3');
subplot (2,2,4)
try
    plot (0:10:LICKS(end),rt);
    axis ([0 LICKS(end) 0 max(rt)]);
catch
end
title ('Licking rate');

%%%% percent events with aLicks
n=0;
for i=1:length(EVENTS{2})
    if sum(LICKS>EVENTS{2}(i)& LICKS<(EVENTS{2}(i)+2))>0
        n=n+1;
    end
end
rez.csp_p=n/length(EVENTS{2});

n=0;
for i=1:length(EVENTS{3})
    if sum(LICKS>EVENTS{3}(i)& LICKS<(EVENTS{3}(i)+2))>0
        n=n+1;
    end
end
rez.csm_p=n/length(EVENTS{3});

%%%%% rate of aLicks
r=zeros(length(EVENTS{2}),1);
for i=1:length(EVENTS{2})
    a1=sum(LICKS>EVENTS{2}(i)-2& LICKS<(EVENTS{2}(i)));
    a2=sum(LICKS>EVENTS{2}(i)& LICKS<(EVENTS{2}(i)+2));
    r(i)=a2/a1;
end
r=r(~isnan(r));
r(r==inf)=10;
rez.csp_r=r;

r=zeros(length(EVENTS{3}),1);
for i=1:length(EVENTS{3})
    a1=sum(LICKS>EVENTS{3}(i)-2& LICKS<(EVENTS{3}(i)));
    a2=sum(LICKS>EVENTS{3}(i)& LICKS<(EVENTS{3}(i)+2));
    r(i)=a2/a1;
end
r=r(~isnan(r));
r(r==inf)=10;
rez.csm_r=r;