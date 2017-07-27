function status=BHV_train5(ct,dly,RWD,F1,F2,F3,N1,N2,N4)
%ct conditioning tone
%dly delay between ct and reward
%RWD amount of reward
%F1,F2,F3 tone frequencies
%N1 minimum ISI
%N2 maximum ISI
%N4 tone duration

% function status=BHV_train5(RWD,F1,F2,F3,N1,N2,N4)

%try,
global LICKS;
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

tone3=sin(F3*[0:(1/Fs):N4]);
tone3(1:length(taps))=tone3(1:length(taps)).*taps;
tone3(length(tone3)-length(tape)+1:end)=tone3(length(tone3)-length(tape)+1:end).*tape;
tone{3}=[tone3' tone3']; % stereo

rand_tones_inx=floor(4*(rand(NTRIALS,1)));
test=1;
while test
    rand_tones_inx(1:100)'
    t=questdlg ('OK ?');
    if t(1)=='Y'
        test=0;
    else
        rand_tones_inx=floor(4*(rand(NTRIALS,1)));
        test=1;
    end

end
rand_times_inx=N1+floor((N2-N1)*(rand(NTRIALS,1)));
rand_tones_inx=ones(NTRIALS,1)*ct;
%try,
    
% initiate the das1200 card 
%CARD('Initialize');
% CARD('lever_up');

%CARD('lever_bit',24);
inx=1;
nrw(1:5)=0;
prog=0.5;
while  inx<NTRIALS,
    nrw(rand_tones_inx(inx)+1)=nrw(rand_tones_inx(inx)+1)+1;
    disp (sprintf ('Tone1=%d \t Tone2=%d \t Tone3=%d \t Catch=%d',nrw(2),nrw(3),nrw(4),nrw(1)));
    if rand_tones_inx(inx)==0
        CARD('send_dio',24);
        putdata(AOS,tone{ct});
        start(AOS);
    elseif rand_tones_inx(inx)==ct
        CARD('send_dio',20+ct);
%         putdata(AOS,tone{ct});
%         start(AOS);
        pause(dly);
        CARD('send_dio',20);
        CARD('send_stimulation');
    else
        CARD('send_dio',20+rand_tones_inx(inx));
        putdata(AOS,tone{rand_tones_inx(inx)});
        start(AOS);
    end
    pause(rand_times_inx(inx));
    inx=inx+1;
end

CARD('lever_down');
CARD('destroy');
close (f);

return;

