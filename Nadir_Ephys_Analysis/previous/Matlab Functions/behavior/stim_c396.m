N1=3;
N2=9;
NTRIALS=1000;

rand_times_inx=N1+floor((N2-N1)*(rand(NTRIALS,1)));
for i=1:100
    disp (i);
    CARD('send_dio',20);
    CARD('send_stimulation');
    pause (rand_times_inx(i));
end


% 
% rand_inxs=floor(2*rand(NTRIALS,1));
% go=1;
% global LICKS
% i=1;
% while go==1
%     disp (['Trial ' num2str(i) ', Int=' num2str(rand_times_inx(i))]);
%     LICKS=0;
%     temp=clock;
%     while LICKS==0
%         if etime(clock, temp)>=rand_times_inx(i)
%             break;
%         end
%         pause (0.001);
%     end
%     if LICKS>0
%         LICKS=0;
%         CARD('send_dio',20);
%         CARD('send_stimulation');
%         pause(dly);
%         rand_times_inx(i)=rand_times_inx(i)-etime(clock, temp);
%         if rand_times_inx(i)<=0
%             CARD('send_dio',21);
%             CARD('reward',0.1);
%             pause(5);
%             i=i+1;
%         end
%     else
%         CARD('send_dio',21);
%         CARD('reward',0.1);
%         pause(5);
%         i=i+1;
%     end
% end