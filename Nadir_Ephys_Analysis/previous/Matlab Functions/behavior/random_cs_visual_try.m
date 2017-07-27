function random_cs_visual_try(REWARD_TIME,P,F1,F2,ITI_BASE,ITI_MAX,swind)

if nargin<1,
    REWARD_TIME=0.1;
end
% pitch frequency
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
    swind=0;
end

scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[1 1 1]);

% initiate the das1200 card 
%CARD('Initialize');

if swind
%    CARD('send_dio',23);
end
pause(ITI_BASE);
for i=1:1000
    p=rand;
    if p<P, % cs plus
        set(gcf,'color',F1);
%        CARD('send_dio',21);
    else
        set(gcf,'color',F2);
%        CARD('send_dio',22);
    end
    pause(1.5);
    set(gcf,'color',[1 1 1]);
    if p<P, % reward
%        CARD('send_dio',20);
%        CARD('reward',REWARD_TIME);
    end
    iti=round(ITI_BASE+rand*ITI_MAX)
    pause(iti);
end

%CARD('destroy');

return;

