function reward(REWARD_TIME)

if nargin<1,
    REWARD_TIME=0.3;
end

CARD('Initialize');

ind=1;
while (ind)
    try,

    pause;    
    CARD('send_dio',20);
    CARD('reward',REWARD_TIME);
    
    catch,
        ind=0;
    end
end

CARD('destroy');

return;

