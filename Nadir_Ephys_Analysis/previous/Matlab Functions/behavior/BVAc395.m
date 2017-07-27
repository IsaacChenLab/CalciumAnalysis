function status=BHV_train5(ct,dly,RWD,F1,F2,N1,N2,N4)
%ct conditioning tone
%dly delay between ct and reward
%RWD amount of reward
%F1,F2,F3 tone frequencies
%N1 minimum ISI
%N2 maximum ISI
%N4 tone duration

% function status=BHV_train5(RWD,F1,F2,F3,N1,N2,N4)

%try,

if nargin<7,
    error('input args');
end
NTRIALS=1e3;

close all;

% prepare the sounds
AOS = analogoutput('winsound');
addchannel(AOS,1:2);
Fs=8e3;
xs=-pi:(pi/160):0;
taps=(cos(xs)+1)/2;
xs=0:(pi/160):pi;
tape=(cos(xs)+1)/2;

tone1=sin(F1*[0:(1/Fs):N4]);
tone1(1:length(taps))=tone1(1:length(taps)).*taps;
tone1(length(tone1)-length(tape)+1:end)=tone1(length(tone1)-length(tape)+1:end).*tape;
tone{1}=[tone1' tone1']; % stereo

tone2=sin(F2*[0:(1/Fs):N4]);
tone2(1:length(taps))=tone2(1:length(taps)).*taps;
tone2(length(tone2)-length(tape)+1:end)=tone2(length(tone2)-length(tape)+1:end).*tape;
tone{2}=[tone2' tone2']; % stereo

rand_tones_inx=floor(2*(rand(NTRIALS,1)))+1;
test=1;
%rand_tones_inx([10 30 50 70 90 110 130 150 170 190])=5;
while test
    rand_tones_inx(1:100)'
    t=questdlg ('OK ?');
    if t(1)=='Y'
        test=0;
    else
        rand_tones_inx=floor(2*(rand(NTRIALS,1)))+1;
        test=1;
        %rand_tones_inx([10 30 50 70 90 110 130 150 170 190])=5;
    end

end
rand_times_inx=N1+floor((N2-N1)*(rand(NTRIALS,1)));

%try,
    
% initiate the das1200 card 
CARD('Initialize');
CARD('lever_up');

CARD('lever_bit',24);
inx=1;
while  inx<NTRIALS,
    disp (rand_tones_inx(inx));
    if rand_tones_inx(inx)==ct 
        CARD('send_dio',20+ct);
        for i=1:5
        putdata(AOS,tone{ct});
        start(AOS);
        pause(3)
        end
%         LICKS=0; %new
        for i=1:4
        putdata(AOS,tone{ct});
        start(AOS);
        pause(1)
        end
%         if LICKS~=0; %new
        CARD('send_dio',20);
        CARD('reward',RWD);
%         end
    else
        CARD('send_dio',20+rand_tones_inx(inx));
        for i=1:6
        putdata(AOS,tone{rand_tones_inx(inx)});
        start(AOS);
        pause(3)
        end
    end        
    pause(rand_times_inx(inx));
    inx=inx+1;

    
end

% catch,
%     disp(lasterr);
% end

CARD('lever_down');
CARD('destroy');

return;

