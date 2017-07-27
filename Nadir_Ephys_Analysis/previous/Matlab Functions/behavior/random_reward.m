function random_reward(REWARD_TIME,P,F1,ITI_BASE,ITI_MAX)

if nargin<1,
    REWARD_TIME=0.1;
end
% pitch frequency
if nargin<3,
    F1=[5e3];
end
if nargin<2
    P=0.8;
end
if nargin<4
    ITI_BASE=10;
end
if nargin<5
    ITI_NAX=10;
end

lngh=0.2;
N1=lngh;
AO = analogoutput('winsound');
addchannel(AO,1:2);
Fs=8e3;
tone=sin(F1*[0:(1/Fs):N1]);
tone=[tone' tone']; % stereo

% initiate the das1200 card 
CARD('Initialize');
    
for i=1:1000
    p=rand;
    putdata(AO,tone);
    CARD('send_dio',21);
    start(AO);
    while strcmp(AO.Running,'On')
        ;
    end
    pause(1);
    if p<P, % reward
        CARD('send_dio',20);
        CARD('reward',REWARD_TIME);
    end
    iti=round(ITI_BASE+rand*ITI_MAX)
    pause(iti);
end

CARD('destroy');

return;

