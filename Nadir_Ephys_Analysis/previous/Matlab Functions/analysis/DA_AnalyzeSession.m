function rez=DA_AnalyzeSession(fl, md, eva, calb)


if nargin<4
    calb=[];
end
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

for i=1:5
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

clear EV EVENTS F Ff H blm fl blp ev evv f ffl fx i j legend lm lp mx n n1 n2 rt x1 x2
[sel k]=listdlg('PromptString','Choose a method:','SelectionMode','single','ListString',{'PCA'; 'Value at +0.7V';...
    'Average +0.6V - +0.8V'; 'Normal distribution average'});
load (ffile);
s=size(SWEEP_DATA);

%%% to see:
% shape of licks outside stimulus/reward
% randomized DA for licks around events

LKS_rnd=sparse(zeros(10000,1));
LKS_rnd(round(1000*EVS{1}))=1;
for i=2:4
    for j=1:length(EVS{i})
        p1=round(1000*max(0.001,EVS{i}(j)-1));
        p2=round(1000*(EVS{i}(j)+5));
        LKS_rnd(p1:p2)=0;
    end
end

% CTR=mean(SWEEP_DATA(985,:))*SWEEP_DATA(4200,:)/mean(SWEEP_DATA(4200,:));
tr=repmat(TRI(1:2000,2),1,size(SWEEP_DATA,2));
if ~isempty (calb)
    CTR=SWEEP_DATA(1:2000,:)./(calb*tr);
else
    CTR=SWEEP_DATA(1:2000,:)./tr;
end
CTR=mean(CTR([900:1300 1500:1950],:),1);
LKS_rnd=find(LKS_rnd)/1000;

% c1=(SWEEP_DATA(4020,i)-SWEEP_DATA(4000,i))/(TRI(4005,1)-TRI(5000,1));
if calib
    cv=abs(DA{1}(4000,10)-DA{1}(2003,10));
    for i=1:s(2)
        SWEEP_DATA(:,i)=(SWEEP_DATA(:,i)*cv)/abs(SWEEP_DATA(4000,i)-SWEEP_DATA(2003,i));
        SWEEP_DATA(:,i)=SWEEP_DATA(:,i)-SWEEP_DATA(4000,i)+DA{1}(4000,10);
    end
    switch sel
        case 1
            
        case 2
            dat=SWEEP_DATA(600,:);
            datf=[];
            for j=1:length(CON)
                datf(j)=mean(DA{j}(600,:));
            end
        case 3
            dat=mean(SWEEP_DATA(800:900,:),1);
            datf=[];
            for j=1:length(CON)
                datf(j)=mean(mean(DA{j}(800:900,:)));
            end
        case 4
            sh=pdf('norm',1:2000,850,15);
            s=size(SWEEP_DATA);
            for i=1:s(2)
                dat(i)=dot(SWEEP_DATA(1:2000,i),sh);
            end
            datf=[];
            for j=1:length(CON)
                datf(j)=mean(dot(DA{j}(1:2000,:),sh,1));
            end
    end
    f=find(CON<300);
    ff=fit(CON(f)',datf(f)','poly1');
    dat=feval(ff,dat);
else
    switch sel
        case 1
            %do PCA

        case 2
            dat=SWEEP_DATA(655,:);
            if ~isempty(calb)
                dat=dat./calb;
            end
        case 3
            dat=mean(SWEEP_DATA(600:700,:),1);
        case 4
            sh=pdf('norm',1:2000,655,15);
            s=size(SWEEP_DATA);
            for i=1:s(2)
                dat(i)=dot(SWEEP_DATA(1:2000,i),sh);
            end
    end
end

q=questdlg('Filter data?','Question','Yes','No','No');
if strcmp(q,'Yes')
    scrsz = get(0,'ScreenSize');
    f=figure('Position',[1*round(scrsz(3)/10) 3*round(scrsz(4)/10) 6*round(scrsz(3)/10) 5*round(scrsz(4)/10)]);
    rep=1;
    pss=0.1;
    while rep
        subplot (2,1,1)
        hold off
        pas=pss/5;
        hd = design(fdesign.highpass ('N,Fp,Ap',3,pas,0.1));
        dd=filter(hd,dat);
        plot (SWEEP_TIMES,dd,'r');
        hold on;
        p1=find (SWEEP_TIMES>=MN);
        p2=find(SWEEP_TIMES>=MX); if isempty(p2) p2=length(SWEEP_TIMES); end
        plot (SWEEP_TIMES,dat-mean(dat(p1(1):p2(1))),'k');
        line([MN MN],[min(dd) max(dd)],'Color','b');
        line([MX MX],[min(dd) max(dd)],'Color','b');
        subplot(2,1,2);
        hold off
        p1=round(length(SWEEP_TIMES)/2);
        plot (SWEEP_TIMES(p1:p1+300),dd(p1:p1+300),'r');
        hold on
        plot (SWEEP_TIMES(p1:p1+300),dat(p1:p1+300)-mean(dat(p1:p1+300)),'k');
        title (['High Pass > ' num2str(pss) 'Hz']);
        [~,~,b]=ginput(1); if isempty (b) b=0; end
        switch b
            case 1 %left
                if pss>=0.02
                    pss=pss-0.01;
                end
            case 3 %right
                pss=pss+0.01;
            otherwise
                rep=0;
                dat=dd;
        end
    end
    close (f);
end

figure;
hold on;
plot (SWEEP_TIMES,dat,'k');
for i=1:length(EVS{1})
    f=find (SWEEP_TIMES>=EVS{1}(i));
    plot (EVS{1}(i),dat(f(1)),'bx');
end
for i=1:length(EVS{2})
    f=find (SWEEP_TIMES>=EVS{2}(i));
    plot (EVS{2}(i),dat(f(1)),'ko');
end
for i=1:length(EVS{3})
    f=find (SWEEP_TIMES>=EVS{3}(i));
    plot (EVS{3}(i),dat(f(1)),'g^');
end
for i=1:length(EVS{4})
    f=find (SWEEP_TIMES>=EVS{4}(i));
    plot (EVS{4}(i),dat(f(1)),'r^');
end

pause(0);
figure;

tms=EVS{1};
df=diff(tms);
f=find(df>=2);
LLKS=intersect(tms(f),tms(f+1));
[shp se]=DA_Avg (dat,SWEEP_TIMES,LLKS,[-5 10]);
rez.DA_LKS=shp;
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,LLKS,[-5 10]);
shpc=shpc-mean(shpc(20:50));
try
    plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
    axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 -0.02 0.05]);
catch end
title ('LKS');

figure;
[shp se]=DA_Avg (dat,SWEEP_TIMES,EVS{1},[-5 10]);
rez.DA_LKS=shp;
subplot (3,3,1);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,EVS{1},[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 -0.02 0.05]);
title ('LKS');


try
[shp se]=DA_Avg (dat,SWEEP_TIMES,EVS{2},[-5 10]);
subplot (3,3,2);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,EVS{2},[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('RWD');
[shpx sex]=rand_lks (EVS{2}, EVS{1}, LKS_rnd, SWEEP_TIMES, dat, [-5 10]);
shpx=shpx-mean(shpx(20:50));
plot (-4.9:0.1:10,shpx,'g'); plot (-4.9:0.1:10,shpx-sex,'g--'); plot (-4.9:0.1:10,shpx+sex,'g--');
rez.DA_RWD=[shp; se; shpc; shpx];
catch
end
try

[shp se]=DA_Avg (dat,SWEEP_TIMES,EVS{3},[-5 10]);
rez.DA_CSP=[shp; se];
subplot (3,3,3);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,EVS{3},[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('CSP');

[shp se]=rand_lks (EVS{3}, EVS{1}, LKS_rnd, SWEEP_TIMES, dat, [-5 10]);
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp,'g'); plot (-4.9:0.1:10,shp-se,'g--'); plot (-4.9:0.1:10,shp+se,'g--');
catch
end
try
[shp se]=DA_Avg (dat,SWEEP_TIMES,EVS{4},[-5 10]);
rez.DA_CSM=[shp; se];
subplot (3,3,4);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,EVS{4},[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('CSM');

[shp se]=rand_lks (EVS{4}, EVS{1}, LKS_rnd, SWEEP_TIMES, dat, [-5 10]);
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp,'g'); plot (-4.9:0.1:10,shp-se,'g--'); plot (-4.9:0.1:10,shp+se,'g--');

catch
end
try
evp=[]; evn=[];
for i=1:length(EVS{3})
    fl=EVS{2}(EVS{2}>EVS{3}(i));
    fl=fl(fl<EVS{3}(i)+4);
    if isempty(fl)
        evn(end+1)=EVS{3}(i);
    else
        evp(end+1)=EVS{3}(i);
    end
end
catch
end
try
[shp se]=DA_Avg (dat,SWEEP_TIMES,evp,[-5 10]);
rez.DA_CSPL=[shp; se];
subplot (3,3,5);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,evp,[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('CSP rewarded');

[shp se]=rand_lks (evp, EVS{1}, LKS_rnd, SWEEP_TIMES, dat, [-5 10]);
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp,'g'); plot (-4.9:0.1:10,shp-se,'g--'); plot (-4.9:0.1:10,shp+se,'g--');
catch
end
try
[shp se]=DA_Avg (dat,SWEEP_TIMES,evn,[-5 10]);
rez.DA_CSPN=[shp; se];
subplot (3,3,6);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,evn,[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('CSP unrewarded');

[shp se]=rand_lks (evn, EVS{1}, LKS_rnd, SWEEP_TIMES, dat, [-5 10]);
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp,'g'); plot (-4.9:0.1:10,shp-se,'g--'); plot (-4.9:0.1:10,shp+se,'g--');
catch
end
try
evp=[]; evn=[];
for i=1:length(EVS{4})
    fl=EVS{1}(EVS{1}>EVS{4}(i)+1);
    fl=fl(fl<EVS{4}(i)+2);
    if isempty(fl)
        evn(end+1)=EVS{4}(i);
    else
        evp(end+1)=EVS{4}(i);
    end
end
catch
end
try
[shp se]=DA_Avg (dat,SWEEP_TIMES,evp,[-5 10]);
rez.DA_CSML=[shp; se];
subplot (3,3,7);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,evp,[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('CSM with LKS');

[shp se]=rand_lks (evp, EVS{1}, LKS_rnd, SWEEP_TIMES, dat, [-5 10]);
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp,'g'); plot (-4.9:0.1:10,shp-se,'g--'); plot (-4.9:0.1:10,shp+se,'g--');
catch
end
try
[shp se]=DA_Avg (dat,SWEEP_TIMES,evn,[-5 10]);
rez.DA_CSMN=[shp; se];
subplot (3,3,8);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
[shpc sec]=DA_Avg (CTR,SWEEP_TIMES,evn,[-5 10]);
shpc=shpc-mean(shpc(20:50));
plot (-4.9:0.1:10,shpc,'r'); plot (-4.9:0.1:10,shpc-sec,'r--'); plot (-4.9:0.1:10,shpc+sec,'r--');
axis ([-5 10 min(shp) max(shp)]);
% axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('CSM w/o LKS');

[shp se]=rand_lks (evn, EVS{1}, LKS_rnd, SWEEP_TIMES, dat, [-5 10]);
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp,'g'); plot (-4.9:0.1:10,shp-se,'g--'); plot (-4.9:0.1:10,shp+se,'g--');
catch
end
try

[shp se]=DA_Avg (dat,SWEEP_TIMES,LKS_rnd,[-5 10]);
rez.DA_CSPL=[shp; se];
subplot (3,3,9);
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('LKS rnd');
catch
end

try

    
drr=EVS{5};
f=zeros(length(drr),1);
for i=1:length(drr)
    if isempty(find(EVS{3}>drr(i)-3 & EVS{3}<drr(i)+3)) && isempty(find(EVS{4}>drr(i)-3 & EVS{4}<drr(i)+3))
        f(i)=1;
    end
end
drr=drr(find(f));
[shp se]=DA_Avg (dat,SWEEP_TIMES,drr,[-5 10]);
rez.DA_CSPL=[shp; se];
figure
hold on
shp=shp-mean(shp(20:50));
plot (-4.9:0.1:10,shp); plot (-4.9:0.1:10,shp-se,'b--'); plot (-4.9:0.1:10,shp+se,'b--');
axis ([-5 10 min(shp)-max(se) max(shp)+max(se)]);
% axis ([-5 10 -0.02 0.05]);
title('LASER');
catch
end


figure;
try
[shp ~]=DA_Avg (SWEEP_DATA(1:2000,:),SWEEP_TIMES,EVS{5},[-5 10]);
catch
    [shp ~]=DA_Avg (SWEEP_DATA(1:2000,:),SWEEP_TIMES,EVS{1},[-5 10]);
end
av=mean(shp(:,20:50),2);
for i=1:150
    shp(:,i)=shp(:,i)-av;
end
LKS_shp=shp;
[shp ~]=DA_Avg (SWEEP_DATA(1:2000,:),SWEEP_TIMES,EVS{2},[-5 10]);
av=mean(shp(:,20:50),2);
for i=1:150
    shp(:,i)=shp(:,i)-av;
end
RWD_shp=shp;
[shp ~]=DA_Avg (SWEEP_DATA(1:2000,:),SWEEP_TIMES,EVS{3},[-5 10]);
av=mean(shp(:,20:50),2);
for i=1:150
    shp(:,i)=shp(:,i)-av;
end
CSP_shp=shp;
[shp ~]=DA_Avg (SWEEP_DATA(1:2000,:),SWEEP_TIMES,EVS{4},[-5 10]);
av=mean(shp(:,20:50),2);
for i=1:150
    shp(:,i)=shp(:,i)-av;
end
CSM_shp=shp;
mnl=min(LKS_shp(:));
mnr=min(RWD_shp(:));
mnp=min(CSP_shp(:));
mnm=min(CSM_shp(:));
mn=min([mnl mnr mnp mnm]);
mnl=max(LKS_shp(:));
mnr=max(RWD_shp(:));
mnp=max(CSP_shp(:));
mnm=max(CSM_shp(:));
mx=max([mnl mnr mnp mnm]);

subplot (2,2,1);
% imagesc (1:150,1:2000,LKS_shp,[mn mx]);
imagesc (1:150,1:2000,LKS_shp);
hold on;
plot ([49 49],[1 2000],'w--', 'LineWidth',3);
title ('LASER');


subplot (2,2,2);
% imagesc (1:150,1:2000,RWD_shp,[mn mx]);
imagesc (1:150,1:2000,RWD_shp);
hold on
plot (2000-2000*ev1(1:150)/max([ev3; ev2; ev1]),'k');
plot ([49 49],[1 2000],'w--', 'LineWidth',3);
title ('RWD');


subplot (2,2,3);
% imagesc (1:150,1:2000,CSP_shp,[mn mx]);
imagesc (1:150,1:2000,CSP_shp);
hold on
plot (2000-2000*ev2(1:150)/max([ev3; ev2; ev1]),'k')
plot ([49 49],[1 2000],'w--', 'LineWidth',3);
plot ([59 59],[1 2000],'w--', 'LineWidth',3);
title ('CSP');


subplot (2,2,4);
% imagesc (1:150,1:2000,CSM_shp,[mn mx]);
imagesc (1:150,1:2000,CSM_shp);
hold on
plot (2000-2000*ev3(1:150)/max([ev3; ev2; ev1]),'k')
plot ([49 49],[1 2000],'w--', 'LineWidth',3);
plot ([59 59],[1 2000],'w--', 'LineWidth',3);
title ('CSM');

CSP_mat=zeros(length(EVS{3}),150);
for i=1:length(EVS{3})
    p=find(SWEEP_TIMES>EVS{3}(i));
    p=p(1);
    CSP_mat(i,:)=dat(p-49:p+100);
end
rez.DA_CSPmat=CSP_mat;
CSM_mat=zeros(length(EVS{4}),150);
for i=1:length(EVS{4})
    p=find(SWEEP_TIMES>EVS{4}(i));
    p=p(1);
    try
        CSM_mat(i,:)=dat(p-49:p+100);
    catch
    end
end
rez.DA_CSMmat=CSM_mat;

figure;
for i=1:length(EVS{3})
    CSP_mat(i,:)=CSP_mat(i,:)-mean(CSP_mat(i,1:49));
end
for i=1:length(EVS{4})
    CSM_mat(i,:)=CSM_mat(i,:)-mean(CSM_mat(i,1:49));
end
% CSP_mat=MC_smooth(MC_smooth(CSP_mat,0.8)',0.8)';
v=CSP_mat(:);
av=mean(v);
st=std(v);
subplot (2,1,1)
imagesc(-4.9:0.1:10,1:length(EVS{3}),CSP_mat,[av-st*4 av+st*4]);
title ('CSP');

% CSM_mat=MC_smooth(MC_smooth(CSM_mat,0.8)',0.8)';
v=CSM_mat(:);
av=mean(v);
st=std(v);
subplot (2,1,2)
imagesc(-49:100,1:length(EVS{4}),CSM_mat,[av-st*4 av+st*4]);
title ('CSM');

end

function [shp se]=rand_lks (evs, lks, lr, st, sd, wnd)

s=size(sd);
shp=zeros(s(1),(wnd(2)-wnd(1))*10);
n=0;
for i=1:length(evs)
    p=round(10*(lks(lks>evs(i)+wnd(1) & lks<evs(i)+wnd(2))-(evs(i)+wnd(1))));
    if ~isempty(p)
        tmp=zeros(s(1),(wnd(2)-wnd(1))*10);
        nn=0;
        for j=1:length(p)
            rep=1;
            while rep
                rep=0;
                xx=lr(1+floor(rand*(length(lr))));
                xxp=find(st>xx);
                xxp=xxp(1);
                try
                    xxd=sd(xxp-(wnd(2)-wnd(1))*10:xxp+(wnd(2)-wnd(1))*10);
                catch
                    rep=1;
                end
            end
            tmp=tmp+xxd((wnd(2)-wnd(1))*10-p(j)+1:(wnd(2)-wnd(1))*10-p(j)+1+(wnd(2)-wnd(1))*10-1);
            nn=nn+1;
        end
        tmp=tmp/nn;
        rep=1;
        while rep
            xp=find(tmp==0);
            if isempty(xp)
                rep=0;
            else
                if xp(1)==1
                    xd=diff(xp);
                    xxp=find (xd>1);
                    if isempty(xxp)
                        tmp(xp)=tmp(xp(end)+1);
                    else
                        tmp(1:xxp(1)+1)=tmp(xxp(1)+2);
                    end
                elseif xp(end)==150
                    tmp(xp(1):end)=tmp(xp(1)-1);
                else
                    xd=diff(xp);
                    xxp=find (xd>1);
                    if isempty(xxp)
                        tmp(xp)=tmp(xp(1)-1);
                    else
                        tmp(xp(1):xp(1+xxp(1)))=tmp(xp(1)-1);
                    end
                end
            end
        end
    else
        tmp=zeros(s(1),(wnd(2)-wnd(1))*10);
    end
    shp(i,:)=tmp;
    shp(i,:)=shp(i,:)-mean(shp(i,1:0-wnd(1)*10));
    n=n+1;
end
se=std(shp,0,1)/sqrt(n);
shp=mean(shp,1);
end