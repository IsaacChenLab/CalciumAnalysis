function reward(REWARD_TIME,t1,t2)

if nargin<1,
    REWARD_TIME=0.3;
end
% pitch frequency
idxs=rand(10000,1)*(t2-t1)+t1;

CARD('Initialize');
CARD('lever_up');

% CARD('lever_bit',24);

for i=1:100000
      CARD('send_dio',20);

    pause(idxs(i));
    CARD('reward',REWARD_TIME);
end
CARD('destroy');

return;

