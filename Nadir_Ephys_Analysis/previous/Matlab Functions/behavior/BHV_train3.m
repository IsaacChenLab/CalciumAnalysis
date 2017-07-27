function status=BHV_train3(RWD,ITI,F1,F2,N1,N2,P,RSP)
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

%try,
global EVENTS;
global EVENT_INX;
global SAVE_CLOCK_START;
EVENTS=[];

if nargin<1,
    RWD=0.1;
end
if nargin<2
    ITI=4;
end
if nargin<3
    F1=4000;
end
if nargin<4,
    F2=7000;
end
if nargin<5,
    N1=6;
end
if nargin<6,
    N2=3;
end
if nargin<7,
    P=0.7;
end
if nargin<8
    RSP=0.5;
end

Fs=8e3;
cs_plus=sin(F1*[0:(1/Fs):N1]);
cs_plus=[cs_plus' cs_plus']; % stereo
cs_minus=sin(F2*[0:(1/Fs):N2]);
cs_minus=[cs_minus' cs_minus'];

% initialize the sound card
AO = analogoutput('winsound');
addchannel(AO,1:2);

% initiate the das1200 card 
CARD('Initialize');
CARD('lever_up');

d=1;
run_ind=1;
while run_ind,
    
%     if rand<0.1
%         CARD('reward',RWD);
%     end
    %BHV_signal('start_trial');
    if rand<P, % cs plus trial
        BHV_signal('cs_plus');
        putdata(AO,cs_plus);
        start(AO); % start sound
        CARD('lever_reward_2',RWD); % activate trigger for lever press with automatic reward
        start_time=clock;
        % sound continues for N1 seconds and while lever is not pressed
        while etime(clock,start_time) < N1 & CARD('check_lever'),
            ;
        end
        stop(AO); % stop sound
        % if lever was not pressed, still allow response time
        while CARD('check_lever') & etime(clock,start_time) < N1+RSP,
            ;
        end
        if ~CARD('check_lever'), % lever was pressed
            BHV_signal('success');
        end
        CARD('lever_stop');
    else % cs minus trial
        BHV_signal('cs_minus');
        press=1;
        N2_curr=N2; % for the first time
        while press==1, % loop continues till lever was not pressed during the sound
            putdata(AO,cs_minus);
            start(AO);
            CARD('lever_wait'); % activate trigger for lever press  
            start_time=clock;
            % sound continues for N2 seconds and while lever is not pressed
            while etime(clock,start_time) < N2_curr & CARD('check_lever'),
                ;
            end
            % if lever was not pressed, continue to next trial.
            if etime(clock,start_time) >= N2_curr, % lever was not pressed (we stopped because of time)
                press=0;
                if N2_curr==N2, % first time
                    BHV_signal('success');
                end
            else, % lever was pressed and this is why we stopped
                BHV_signal('press');
            end
            CARD('lever_stop');
            stop(AO);
            N2_curr=round(N2*1.5); % each press lengthen by half of the original
        end
    end
    
    BHV_signal('ITI');
    % inter trial interval
    iti=round(ITI+rand*ITI);
    pause(iti);

    %     if length(ITI)==2
%         IPI=rand(1)*range(ITI) + ITI(1); % generate a random interval
%     else
%         IPI=ITI;
%     end
%     press=1;
%     IPI_curr=IPI;
%     while press==1, % loop continues till lever was not pressed during the ITI
%         CARD('lever_wait'); % activate trigger for lever press  
%         start_time=clock;
%         while etime(clock,start_time) < IPI & CARD('check_lever'),
%             ;
%         end
%         if etime(clock,start_time) >= IPI_curr, % lever was not pressed (we stopped because of time)
%             press=0;
%         else
%             BHV_signal('press');
%         end
%         CARD('lever_stop');
%         %IPI_curr=round(IPI*1.5);
%     end
    d=d+1;
end

CARD('lever_stop');
CARD('lever_down');
CARD('destroy');
BHV_signal('save_events');

% catch,
%     BHV_signal('save_events');
%     CARD('lever_stop');
%     CARD('lever_down');
%     CARD('destroy');
%    
% end

return;

