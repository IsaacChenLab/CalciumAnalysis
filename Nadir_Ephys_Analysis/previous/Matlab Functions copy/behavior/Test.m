function status=BHV_train5(F1,F2,F3,N1,N4)
%ct conditioning tone
%dly delay between ct and reward
%RWD amount of reward
%F1,F2,F3 tone frequencies
%N1 minimum ISI
%N2 maximum ISI
%N4 tone duration

% function status=BHV_train5(RWD,F1,F2,F3,N1,N2,N4)

%try,

NTRIALS=4;

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

tone3=sin(F3*[0:(1/Fs):N4]);
tone3(1:length(taps))=tone3(1:length(taps)).*taps;
tone3(length(tone3)-length(tape)+1:end)=tone3(length(tone3)-length(tape)+1:end).*tape;
tone{3}=[tone3' tone3']; % stereo

rand_tones_inx=[3 1 3 2];
rand_times_inx=[N1 2*N1 3*N1 4*N1];

%try,
    
% initiate the das1200 card 
CARD('Initialize');
CARD('lever_up');

CARD('lever_bit',24);
inx=1;
while  inx<NTRIALS,
    disp (rand_tones_inx(inx));
%     if rand_tones_inx(inx)==ct || rand_tones_inx(inx)==0
%         CARD('send_dio',20+ct);
%         putdata(AOS,tone{ct});
%         start(AOS);
%         pause (dly);
%         CARD('send_dio',20);
%         CARD('reward',RWD);
    if rand_tones_inx(inx)==5
        CARD('send_dio',20+ct);
        putdata(AOS,tone{ct});
        start(AOS);
    else
        CARD('send_dio',20+rand_tones_inx(inx));
        putdata(AOS,tone{rand_tones_inx(inx)});
        start(AOS);
    end        
%     if rand_tones_inx(inx)>0, % tone
%         putdata(AOS,tone{rand_tones_inx(inx)});
%         start(AOS);
%     elseif rand_tones_inx(inx)==0,
%         CARD('reward',RWD); % reward
%     end
    while strcmp(AOS.Running,'On')
        ;
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

