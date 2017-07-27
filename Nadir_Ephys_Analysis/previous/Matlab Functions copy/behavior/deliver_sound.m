function deliver_sound()

flip_flop=0;
lngh=0.01;
which=20;
N1=lngh;
F=[15e3 20e3];
AO = analogoutput('winsound');
addchannel(AO,1:2);
Fs=8e3;

CARD('Initialize');
    F1=F(1+flip_flop);
    tone=sin(F1*[0:(1/Fs):N1]);
    tone=[tone' tone']; % stereo

for i=1:100,
    if rem(i,10)==0
        flip_flop=1-flip_flop;
        F1=F(1+flip_flop);
        tone=sin(F1*[0:(1/Fs):N1]);
        tone=[tone' tone']; % stereo
    end

    putdata(AO,tone);
    CARD('send_dio',which);
    start(AO);
    pause(3);
end

CARD('destroy');

return;
