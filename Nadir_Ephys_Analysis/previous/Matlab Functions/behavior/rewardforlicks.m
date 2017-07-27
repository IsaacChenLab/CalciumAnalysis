function reward(REWARD_TIME,t1,t2)

if nargin<1,
    REWARD_TIME=0.3;
end
ri=floor(REWARD_TIME/0.1);
% pitch frequency
idxs=rand(10000,1)*(t2-t1)+t1;

global LICKS;
% CARD('lever_bit',24);

for i=1:100000

    pause(t1);
    temp=clock;
    LICKS=0;
    pause (0.3);
    while etime(clock, temp)<t2 && LICKS==0
        pause(0.3);
    end
    if LICKS>0
        CARD('send_dio',20);
        for jj=1:ri
            CARD('reward',0.1);
            pause(0.6);
        end
    else
        CARD('send_dio',20);
        for jj=1:ri
            CARD('reward',0.1);
            pause(0.6);
        end
    end
end
CARD('destroy');

return;

