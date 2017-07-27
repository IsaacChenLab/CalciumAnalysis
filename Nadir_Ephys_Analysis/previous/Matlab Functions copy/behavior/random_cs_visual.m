function random_cs(REWARD_TIME,P,F1,F2,ITI_BASE,ITI_MAX,DELAY,SWTIME)

close all;

if nargin<1,
    REWARD_TIME=0.1;
end
% color
if nargin<3,
    F1=[1 0 0];
end
if nargin<4
    F2=[0 1 0];
end
if nargin<2
    P=0.5;
end
if nargin<5
    ITI_BASE=10;
end
if nargin<6
    ITI_MAX=10;
end
if nargin<7
    DELAY=3;
end
if nargin<8
    SWTIME=0;
end

onethird=DELAY/3;
twothirds=DELAY-onethird;

scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[1 1 1]);
    
% initiate the das1200 card 
CARD('Initialize');

try,
    
if SWTIME
    CARD('send_dio',23);
end

pause(ITI_BASE);
ind=1;
d=1;
while (ind)
    try,
    if SWTIME & mod(d,SWTIME)==0
        CARD('send_dio',23);
        F3=F2; F2=F1; F1=F3;
        disp('switching');
    end
    if mod(d,SWTIME)>(SWTIME/2)
        PP=0.9;
    else
        PP=1;
    end
    pp=rand;
    p=rand;
    delay=onethird+(rand*twothirds);
    restdelay=DELAY-delay-REWARD_TIME;
    if p>=P, % cs plus
        set(gcf,'color',F1);
        CARD('send_dio',21);
    else
        set(gcf,'color',F2);
        CARD('send_dio',22);
    end
    pause(delay);
    if p>=P & pp<PP, % reward
        CARD('send_dio',20);
        CARD('reward',REWARD_TIME);
    end
    pause(restdelay);
    set(gcf,'color',[1 1 1]);

    iti=round(ITI_BASE+rand*ITI_MAX);
    pause(iti);
    d=d+1;
    catch,
        ind=0;
    end
end
    
CARD('destroy');
close all;

end

return;

