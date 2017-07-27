function reward(REWARD_TIME,CS_FREQ)

if nargin<1,
    REWARD_TIME=0.3;
end
% pitch frequency
if nargin<2,
    CS_FREQ=1000;
end

CARD('Initialize');

while (1)
    
lngh=0.2;
N1=lngh;
F1=[5e3];
AO = analogoutput('winsound');
addchannel(AO,1:2);
Fs=8e3;
tone=sin(F1*[0:(1/Fs):N1]);
tone=[tone' tone']; % stereo
%putdata(AO,tone);
%CARD('send_dio',21);
%start(AO);
% pause(3);


% initiate the das1200 card 
    CARD('send_dio',20);

    CARD('reward',REWARD_TIME);
    pause;
end

CARD('destroy');

return;

