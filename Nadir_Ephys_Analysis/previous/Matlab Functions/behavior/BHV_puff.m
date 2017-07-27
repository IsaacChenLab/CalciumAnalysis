function status=BHV_train5(RWD,ITI,F1,F2,N1,N2,P,SWTIME)
% function status=BHV_train2(RWD,ITI,F1,F2,N1,N2,P,RSP)
% training paradigm:
% two tone frequencies: F1 , F2
% each trial a tone is selected randomly (with P probability for F1 an 1-P
% for F2), and sounded for N1/N2 seconds.
% if the lever is pressed while F1 is sounded (or within RSP seconds after it ends), then reward is delivered for RWD seconds.
% if the lever is pressed while F2 is sounded, the tone extends for N2 more
% seconds.
% inter trial interval for ITI (can be a range, if so selected random) seconds, if the lever is pressed, ITI
% extends for ITI more seconds.

CARD('Initialize');
CARD('lever_up');

CARD('lever_bit',24);
d=1;
run_ind=1;
while run_ind
    CARD('send_stimulation');
    pause(5);
end

CARD('lever_down');
CARD('destroy');

return;

