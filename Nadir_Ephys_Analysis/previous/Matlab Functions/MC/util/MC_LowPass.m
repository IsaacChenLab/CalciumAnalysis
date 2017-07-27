function fltData=MC_spHighPass(rawData,Fs,hp)
% function fltData=MC_highpass(rawData,Fs)
% filters high pass from 150hz

if nargin<2
    Fs=25e3;
end
if nargin<3,
    hp=300;
end

Nyquist=Fs/2;
Wp=[hp]/Nyquist;
%[B,A]=butter(2,Wp,'high');
%[B,A]=butter(6,Wp,'high');
[B,A]=cheby2(4,35,Wp,'low');
fltData=filter(B,A,rawData);

return;