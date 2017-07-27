
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
    N1=20;
end
if nargin<4,
    N2=20;
end

NTRIALS=1e3;

close all;

% initiate the das1200 card
CARD('Initialize');
CARD('lever_up'); % putvalue(DIO.Line(19),1)

scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);drawnow;

%random
r=rand(NTRIALS,1);
rr=rand(NTRIALS,1);
%inter trial interval
rand_times=N1+floor((N2-N1)*rr);
D1=10;D2=18; %a remettre
rand_times_horizlines=D1+floor((D2-D1)*rr); %a remettre
k=rand_times_horizlines; %a remettre
% times(1:2:10000)=0; times(2:2:10000)=1;
% for i=1:4:9996
%     pos=floor(rand(1,1)*4)+1;
%     times2=[times2; times(i:i+pos-1);times(i+pos-1);times(i+pos:i+3)];
% end

CARD('lever_bit',24);
d=1;
run_ind=1;
while run_ind,
    if  r(d)<0.5, %50% probability;  cs plus trial
        CARD('send_dio',21);
        n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
        t(find(t==0))=10;
        t=t-10;
        %         oriz,for i=1:1000, axis([1 5 t(i) t(i)+10]);, pause(0.01);,end,close
        %         LICKS=0;

        CS_time=k; %CS_time=16 a remettre apres
        RWD_time=4;
        RATIO=2;
        temp=clock;
        oriz
        i=1;
        LICKS=0; %a remettre apres
        while 1==1
            axis([1 5 t(i) t(i)+10]);
            i=i+1;
            if etime(clock, temp)>=CS_time
                break;
            end
            pause(0.01);
        end
        close;
        base=LICKS/CS_time;

        temp=clock;
        oriz1;
        i=1;

        test=0; %a remettre apres
        LICKS=0; % a remettre apres

        n=0;for i=1:0.03:1000, n=n+1;t(n)=i-floor(i/10)*10;,end % m=2
        t(find(t==0))=10;
        t=t-10;
        while 1==1
            axis([1*cos(63) 2*cos(63) t(i)*sin(63) (t(i)+10)*sin(63)]);
            i=i+1;
            if etime(clock, temp)>=RWD_time
                break;
            end
            pause(0.01);
        end
        close;

        if LICKS/RWD_time>base*RATIO %a remettre apres
            CARD('send_dio',20);
            CARD('reward',RWD);
        end % a remettre apres

        %%n=0;for i=1:0.03:1000, n=n+1;t(n)=i-floor(i/10)*10;,end %m=1
        %%t(find(t==0))=10;
        %%t=t-10;
        %%oriz1,for i=1:1000, axis([1*cos(45) 2*cos(45) t(i)*sin(45) (t(i)+10)*sin(45)]);, pause(0.01);,end,close
        %%n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end % m=2
        %%t(find(t==0))=10;
        %%t=t-10;
        %%oriz1,for i=1:632, axis([1 2 t(i)/2 (t(i)+10)/2]);, pause(0.01);,end,close
        %%pause(dly);
        %%CARD('send_dio',20);
        %%CARD('reward',RWD);
        %         n=0;for i=1:0.03:1000, n=n+1;t(n)=i-floor(i/10)*10;,end % m=2
        %         t(find(t==0))=10;
        %         t=t-10;
        %         oriz1,for i=1:167, axis([1*cos(63) 2*cos(63) t(i)*sin(63) (t(i)+10)*sin(63)]);, pause(0.01);,end,close
        %         if LICKS~=0;% licks during sloping lines
        %         %pause(dly);
        %         CARD('send_dio',20);
        %         CARD('reward',RWD);
        %         end
        %no licks, no reward, we continue to next trial
    else % cs minus trial
        CARD('send_dio',22);
        %f1; pause(20); %black screen
        n=0;for i=1:0.0125:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
        t(find(t==0))=10;
        t=t-10;
        %%         vert,for i=1:990, axis([t(i) t(i)+10 0 10]);, pause(0.01);,end,close

        CSminus_time=20;
        vert
        i=1;
        temp=clock;
        while 1==1
            axis([t(i) t(i)+10 0 10]);
            i=i+1;
            if etime(clock, temp)>=CSminus_time
                break;
            end
            pause(0.01);
        end
        close;

    end
    pause(rand_times(d));% inter trial interval - end
    d=d+1;


end
close(f1); % close the initial figure
CARD('lever_stop');
CARD('lever_down'); % putvalue(DIO.Line(19),0)
CARD('destroy');

return;