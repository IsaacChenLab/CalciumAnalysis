
function BV(RWD,dly,N1,N2)
% training paradigm:
% each visual trial is selected randomly with a probability of 50% cs+
% dly delay between cs and reward (s)
% RWD amount of reward (s)
% N1 minimum inter trial interval (s)
% N2 maximum inter trial interval (s)
% january, 26, 2007

if nargin<1,
    RWD=0.15;
end
if nargin<2
    dly=3; % > if we are interested in hippocampus
end
if nargin<3
    N1=20;
end
if nargin<4,
    N2=40;
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
        circle([1,3],2,100000,'w');axis off,
        pause(3);clf;
        circle([1,3],2,100000,'k');axis off; hold on; circle([1,3],.25,100000,'w'); axis off; drawnow;
        pause(dly);clf;
        CARD('send_dio',20);
        CARD('reward',RWD);
    else % cs minus trial  
        CARD('send_dio',22);
        circle([1,3],2,100000,'k');axis off; hold on; circle([1,3],1,100000,'b'); axis off; drawnow;
        pause(3);clf;
    end
    pause(rand_times(d));% inter trial interval - end
    d=d+1;
    
end

close(f1); % close the initial figure
CARD('lever_down'); % putvalue(DIO.Line(19),0)
CARD('destroy');

return;