function status=BHV_train5(RWD,ITI,F1,F2,N1,N2,P,RSP,SWTIME)
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
    F1=[0 1 0];
end
if nargin<4,
    F2=[0 0 1];
end
if nargin<5,
    N1=6;
end
if nargin<6,
    N2=3;
end
if nargin<7,
    P=0.5;
end
if nargin<8
    RSP=0.5;
end
if nargin<9
    SWTIME=0;
end

close all;

% initiate the das1200 card 
CARD('Initialize');
CARD('lever_up');

scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[1 1 1]);
drawnow;

try,

ps=rand(1e4,1);

if SWTIME
    CARD('send_dio',23); CARD('send_dio',23); CARD('send_dio',23);
end

d=1;
run_ind=1;
while run_ind,
    if SWTIME & mod(d,SWTIME)==0
        CARD('send_dio',23); CARD('send_dio',23); CARD('send_dio',23);
        F3=F2; F2=F1; F1=F3;
        disp('switching');
    end

    if ps(d)<P, % cs plus trial
        CARD('send_dio',21);
        set(gcf,'color',F1); drawnow;
        CARD('lever_wait');
        start_time=clock;
        % sound continues for N1 seconds and while lever is not pressed
        while etime(clock,start_time) < N1 & CARD('check_lever'),
            ;
        end
        set(gcf,'color',[1 1 1]); drawnow;
        % if lever was not pressed, still allow response time
        while CARD('check_lever') & etime(clock,start_time) < N1+RSP,
            ;
        end
        if ~CARD('check_lever'), % lever was pressed
            CARD('send_dio',20);
            CARD('reward',RWD);
        end
        CARD('lever_stop');
    else % cs minus trial
        CARD('send_dio',22);
        press=1;
        first=1;
        N2_curr=N2; % for the first time
        while press==1, % loop continues till lever was not pressed during the sound
            set(gcf,'color',F2); drawnow;
            CARD('lever_wait'); % activate trigger for lever press  
            start_time=clock;
            % sound continues for N2 seconds and while lever is not pressed
            while etime(clock,start_time) < N2_curr & CARD('check_lever'),
                ;
            end
            % if lever was not pressed, continue to next trial.
            if etime(clock,start_time) >= N2_curr, % lever was not pressed (we stopped because of time)
                press=0;
                if first, % first time
                        %CARD('send_dio',20);
                end
            else, % lever was pressed and this is why we stopped
                set(gcf,'color',[1 1 1]); drawnow; pause(0.1);
                first=0;
            end
            CARD('lever_stop');
            N2_curr=N2; % each press lengthen by the original
        end
    end
    
    set(gcf,'color',[1 1 1]); drawnow;
    % inter trial interval
    iti=round(ITI+ps(d)*ITI);
    CARD('send_dio',23);
    CARD('lever_wait');
    pause(iti);
    if ~CARD('check_lever')
        %CARD('send_dio',24);
    else
        CARD('lever_stop');
    end
    d=d+1;
end

catch,
    disp(lasterr);
end

close(f1);
CARD('lever_stop');
CARD('lever_down');
CARD('destroy');

return;

