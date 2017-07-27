
function BV(RWD,dly,N1,N2)
% training paradigm:
% each visual trial is selected randomly with a probability of 50% cs+
% dly delay between cs and reward (s)
% RWD amount of reward (s)
% N1 minimum inter trial interval (s)
% N2 maximum inter trial interval (s)
% january, 26, 2007

global LICKS;

if nargin<1,
    RWD=0.15;
end
if nargin<2
    dly=3; % > if we are interested in hippocampus tasks
end
if nargin<3
    N1=3;
end
if nargin<4,
    N2=6;
end

NTRIALS=1e3;

close all;

% initiate the das1200 card 
CARD('Initialize'); 
CARD('lever_up'); % putvalue(DIO.Line(19),1)

scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);
drawnow;

%random
r=rand(NTRIALS,1);
rr=rand(NTRIALS,1);
%inter trial interval
rand_times=N1+floor((N2-N1)*rr);

CARD('lever_bit',24);
d=1;
run_ind=1;
while run_ind,
    if r(d)<0.5, %50% probability;  cs plus trial
        CARD('send_dio',21);
        for i=1:10
        circle([0,0],2,1000,'k');hold on; circle([0,0],.25,1000,'w'); drawnow;
        pause(.3);clf;pause(.2);
        end
        p=rand(1,1);n=round(5*p); %value between 0 and X, here X=5.
        for i=1:n
        circle([0,0],2,1000,'k');hold on; circle([0,0],.25,1000,'w'); drawnow;
        pause(.3);clf;pause(.2);
        end
        %D1=4;D2=10;
        %rand_times_circle=D1+floor((D2-D1)*rr);
        %k=rand_times_circle;
        %for i=1:k
        %circle([1,3],2,1000,'k'); hold on; circle([1,3],.25,1000,'w'); drawnow;
        %pause(.3);clf;pause(.2);
        %end
        LICKS=0;
        for i=1:10
        %circle([0,0],2,1000,'w');
        %pause(.3);clf;pause(.2);
        circle([0,0],2,1000,'k');hold on; circle([0,0],.05,1000,'w'); drawnow;
        pause(.3);clf;pause(.2);
        end
        if LICKS==1;% licks during the big (or very small) circles 
        pause(dly);
        CARD('send_dio',20);
        CARD('reward',RWD);
        end
        %no licks, no reward, we continue to next trial    
    else % cs minus trial  
        CARD('send_dio',22);
        for i=1:10 %i=1:180 for 90s
        circle([0,0],2,1000,'k'); hold on; rectangle2(.5,.1,'w'); drawnow;
        pause(.3);clf;pause(.2);
        end
    end
    pause(rand_times(d));% inter trial interval - end
    d=d+1;
    
end

close(f1); % close the initial figure
CARD('lever_stop');
CARD('lever_down'); % putvalue(DIO.Line(19),0)
CARD('destroy');

return;