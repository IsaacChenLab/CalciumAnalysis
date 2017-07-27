
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
    dly=0;
end
if nargin<3
    N1=0;
end
if nargin<4,
    N2=0;
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
        n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
        t(find(t==0))=10;
        t=t-10;
        oriz,for i=1:1000, axis([1 5 t(i) t(i)+10]);, pause(0.01);,end,close
        LICKS=0;
        %%n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end %m=1
        %%t(find(t==0))=10;
        %%t=t-10;
        %%oriz12,for i=1:632, axis([1 2 t(i)/2 (t(i)+10)/2]);, pause(0.01);,end,close
        %%n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end % m=2
        %%t(find(t==0))=10;
        %%t=t-10;
        %%oriz1,for i=1:632, axis([1 2 t(i)/2 (t(i)+10)/2]);, pause(0.01);,end,close
        %%pause(dly);
        %%CARD('send_dio',20);
        %%CARD('reward',RWD);
        n=0;for i=1:0.03:1000, n=n+1;t(n)=i-floor(i/10)*10;,end % m=2
        t(find(t==0))=10;
        t=t-10;
        oriz1,for i=1:167, axis([1*cos(63) 2*cos(63) t(i)*sin(63) (t(i)+10)*sin(63)]);, pause(0.01);,end,close
        if LICKS==1;% licks during sloping lines
        %pause(dly);
        CARD('send_dio',20);
        CARD('reward',RWD);
        end
        %no licks, no reward, we continue to next trial    
    else % cs minus trial  
        CARD('send_dio',22);
        f1; pause(20);
        %n=0;for i=1:0.0125:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
        %t(find(t==0))=10;
        %t=t-10;
        %vert,for i=1:1430, axis([t(i) t(i)+10 0 10]);, pause(0.01);,end,close
    end
    pause(rand_times(d));% inter trial interval - end
    d=d+1;
    
end

close(f1); % close the initial figure
CARD('lever_stop');
CARD('lever_down'); % putvalue(DIO.Line(19),0)
CARD('destroy');

return;