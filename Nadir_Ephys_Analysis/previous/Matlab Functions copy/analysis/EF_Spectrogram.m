function r = EF_Spectrogram (dat,Fq, freq, base, minwnd)

% dat = raw data
% Fq = frequencies for analisys
% freq = sampling frequency
% base = baseline times
% minwnd = do not go below this value for time window
r.tms=[];
r.pow=[];
r.freq=[];
r(length(freq)).tms=[];
h=waitbar(0);
for i=1:length(Fq)
    wnd=round(max(minwnd*freq,4*freq/Fq(i)));
    [S,F,T,P]=spectrogram(dat,wnd,round(wnd/1.1),[Fq(i)-.1 Fq(i)],freq);
    r(i).tms=T;
    r(i).pow=P(2,:);
    r(i).freq=Fq(i);
    ff=find(r(i).tms>min(base));
    p1=ff(1);
    ff=find(r(i).tms<max(base));
    p2=ff(end);
    r(i).pown=r(i).pow/mean(r(i).pow(p1:p2));
    waitbar(i/length(Fq));
end
close (h);
beep
figure;
hold on;
for i=1:length(Fq)
    imagesc(r(i).tms,r(i).freq,r(i).pown);
end
    