function [F T rez]=EF_MovingCoherence (lfp1, lfp2, wnd, stp, samp)

% calculates coherence of the two lfps (with recording freq samp) in short windows of wnd (sec.)
% length with stp sliding between steps, for frequencies lower than freq
%[Cxy,F] = mscohere(x,y,window,noverlap,nfft,fs)

win=round(samp*0.2);
noo=round(win/1.2);
len=floor(wnd*samp);
stps=round(stp*samp);
[C F]=mscohere(lfp1(1:1+len),lfp2(1:1+len), win, noo,round(samp),samp);
rez=zeros(length(F),1);
T=0;
h=waitbar (0);
for i=1:stps:(length(lfp1)-wnd*samp)
    waitbar(i/(length(lfp1)-wnd*samp));
    [rez(:,end+1) ~]=mscohere(lfp1(i:i+len),lfp2(i:i+len), win, noo,round(samp),samp);
    T(end+1)=T(end)+stp;
end
close(h);
rez=rez(:,2:end);
T=T(1:end-1);