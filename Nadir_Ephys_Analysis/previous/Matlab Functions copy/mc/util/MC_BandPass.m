function fltData=MC_spHighPass(rawData,Fs,bp)
% function fltData=MC_highpass(rawData,Fs)
% filters high pass from 150hz

if nargin<2
    Fs=25e3;
end
if nargin<3,
    bp=[2 150];
end

Nyquist=Fs/2;
Wp=bp/Nyquist;
%[B,A]=butter(2,Wp,'high');
%[B,A]=butter(6,Wp,'high');
[B,A]=cheby2(2,35,Wp,'bandpass');
fltData=filtfilt(B,A,rawData);

return;