function random_cs(REWARD_TIME,P,F1,F2,ITI)

if nargin<1,
    REWARD_TIME=0.1;
end
% pitch frequency
if nargin<3,
    F1=[5e3];
end
if nargin<4
    F2=F1-1e3;
end
if nargin<2
    P=0.5;
end
if nargin<5
    ITI=10;
end

lngh=0.3;
N1=lngh;
AO = analogoutput('winsound');
addchannel(AO,1:2);
Fs=8e3;
tone1=sin(F1*[0:(1/Fs):N1]);
tone1=[tone1' tone1']; % stereo
tone2=sin(F2*[0:(1/Fs):N1]);
tone2=[tone2' tone2']; % stereo

% initiate the das1200 card 
CARD('Initialize');
CARD('lever_up');

CARD('send_dio',23);

while(1),
    
    CARD('lever_cs');
    p=rand;
    if p<P, % cs plus
        putdata(AO,tone1);
        CARD('send_dio',21);
    else
        putdata(AO,tone2);
        CARD('send_dio',22);
    end
    start(AO);
    while strcmp(AO.Running,'On')
        ;
    end
    pause(1.5);
    if p<P, % reward
        CARD('send_dio',20);
        CARD('reward',REWARD_TIME);
    end
    pause(ITI);
end

CARD('lever_down');
CARD('destroy');

return;

