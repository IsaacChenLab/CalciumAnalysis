function rez=DA_AnalyzeBehavior(fl, md, eva)

if nargin<3
    eva=[];
end
if nargin<2
    md='all';
end

tmp=which(fl);
ff=find(tmp=='\');
if exist ([tmp(1:ff(end)) 'CALIBRATION.mat'],'file')
    a=questdlg('Convert data to concentration using CALIBRATION.mat file?','Question','Yes','No','No');
    if strcmp(a,'Yes')
        load ([tmp(1:ff(end)) 'CALIBRATION.mat']);
        calib=1;
    else
        calib=0;
    end
else
    calib=0;
end

rng('shuffle');
rez=[];
load (fl,'EVENTS');
switch lower(md)
    case 'tone'
        EVENTS(3,:)=EVENTS(5,:);
        EVENTS(4,:)=EVENTS(6,:);
    case 'light'
        EVENTS(3,:)=EVENTS(7,:);
        EVENTS(4,:)=EVENTS(8,:);
    case 'all'
    
    otherwise
        error ('not a valid input argument');
end
if isempty(eva)
    rt=Rate(find(EVENTS(1,:))/1000,10)/10;
    f=figure;
    plot (rt);
    hold on;
    [x1 ~]=ginput(1);
    line([x1 x1],[min(rt) max(rt)]);
    [x2 ~]=ginput(1);
    close (f);
    MX=round(10*max(x1,x2));
    MN=round(10*min(x1,x2));
else
    rt=Rate(find(EVENTS(1,:))/1000,10)/10;
    f=figure;
    MX=15*(max(eva)-1);
    MN=15*(min(eva)-1);
    x1=MN;
    x2=MX;
    rectangle ('Position', [MN/10 0 (MX-MN)/10 max(rt)], 'FaceColor','y','EdgeColor','y');
    hold on
    plot (rt);
    ginput(1);
    close (f);
end
ffile=fl;

for i=1:4
    try
        ev=find(EVENTS(i,:))/1000;
        EV=ev(ev<MX);
        EV=EV(EV>MN);
        EVS{i}=EV;
    catch
    end
end

H=0;
F=0;
Ff=0;

for i=1:length(EVS{3})
    fl=find(EVS{1}<EVS{3}(i)+1);
    try if EVS{1}(fl(end))>EVS{3}(i)
        H=H+1;
        end
    catch
    end
end
H=H/length(EVS{3});
if H==1
    H=0.999999;
elseif H==0
    H=0.000001;
end

for i=1:length(EVS{4})
    fl=find(EVS{1}<EVS{4}(i)+1);
    ffl=find(EVS{1}<EVS{4}(i)+2);
    try if EVS{1}(fl(end))>EVS{4}(i)
        F=F+1;
        end
    catch
    end
    try if EVS{1}(ffl(end))>(EVS{4}(i)+1)
        Ff=Ff+1;
        end
    catch
    end
end
F=F/length(EVS{4});
if F==1
    F=0.999999;
elseif F==0
    F=0.000001;
end
Ff=Ff/length(EVS{4});
if Ff==1
    Ff=0.999999;
elseif Ff==0
    Ff=0.000001;
end

for i=2:4
    try
        eval (['ev' num2str(i-1) '=MC_smooth(XCorr (EVS{' num2str(i) '},EVS{1},[-5 10],.1)'',0.5);']);
    catch
    end
end

mx=0;
for i=1:3
    try
        mx=max(mx,max(eval(['ev' num2str(i)])));
    catch
    end
end
figure

legend={'REWARD'; 'CSP'; 'CSM'};
for i=1:3
    try
        subplot (2,2,i);
        plot (-5:.1:10,eval(['ev' num2str(i)]));
        axis ([-5 10 0 mx]);
        title (legend{i});
    catch
    end
end

if exist('ev3','var')
    subplot (2,2,4)
else
    subplot(2,2,3:4)
end

fx=find(EVENTS(1,:))/1000;
plot (0:10:fx(end),rt);
hold on
axis ([0 EVS{1}(end) 0 max(rt)]);
line([x1*10 x1*10],[min(rt) max(rt)],'Color','r');
line([x2*10 x2*10],[min(rt) max(rt)],'Color','r');
title ('Licking rate');

evv=[];
for i=1:length(EVS{2})
    if ~isempty(find(EVS{1}>EVS{2}(i) & EVS{1}<EVS{2}(i)+1))
        evv(end+1)=EVS{2}(i);
    end
end
EVS{2}=evv;

for i=3:4
    try
        ev=eval(['EVS{' num2str(i) '}']);
        n=0;
        for j=1:length(ev)
            if ~isempty (find(EVS{2}>ev(j) & EVS{2}<ev(j)+5))
                n=n+1;
            end
        end
        eval(['rez.ev' num2str(i) '=n/length(ev);']);
        eval (['n' num2str(i-2) '=n;']);
    catch
    end
end

rez.cor=n1/length(EVS{3});
lp=sum(ev2(51:62));
blp=sum(ev2(41:50));
lm=sum(ev3(51:62));
blm=sum(ev3(41:50));
rez.percent_lp=lp/blp;
rez.percent_lm=lm/blm;
% lp=max(0,lp-blp);
% lm=max(0,lm-blm);
rez.dis=(lp-lm)/(lp+lm);
rez.dprime_an=norminv(H,0,1)-norminv(F,0,1);
rez.dprime_po=norminv(max(min(rez.cor,0.999999),0.000001),0,1)-norminv(Ff,0,1);
rez.false=Ff;