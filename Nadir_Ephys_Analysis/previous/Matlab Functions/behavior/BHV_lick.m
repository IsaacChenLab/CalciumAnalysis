function status=BHV_lick(RWD)
if nargin==0
    RWD=0.12;
end


CARD('Initialize');
CARD('lever_up');

CARD('lever_bit',24);
run_ind=1;
while run_ind,

    waitforbuttonpress;
    CARD('reward',RWD);

    % inter trial interval
    pause(1);
end


CARD('lever_down');
CARD('destroy');

return;

