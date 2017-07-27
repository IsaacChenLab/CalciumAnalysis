function status=BHV_train2(REWARD_TIME,ITI,CS_FREQ)

if nargin<1,
    REWARD_TIME=0.5;
end
if nargin<2
    ITI=2;
end
% pitch frequency
if nargin<3,
    CS_FREQ=1000;
end

% initiate the das1200 card 
CARD('Initialize');


CARD('lever_up');

run_ind=1;
while run_ind,
    
    CARD('lever_sound_reward',REWARD_TIME,CS_FREQ);
    pause(ITI);
    
end

CARD('lever_down');

CARD('destroy');

return;

