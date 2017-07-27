function EF_PlotRaw (fl, evs, wnd, flsub, patt)

if nargin < 5
    patt=[];
end
if nargin < 4
    flsub=[];
end
load (fl)
r1=raw;
if ~isempty(flsub)
    load (flsub);
    r0=raw;
    clear raw;
else
    r0=r1;
    r0(:)=0;
end
nn=length(evs);
tmp=abs(r1-r0);
mn=mean(tmp);
sd=std(tmp);
[~, spk]=findpeaks(double(tmp),'MINPEAKHEIGHT',mn+4*sd);
clear tmp;
for i=1:nn
    eval(['[m' num2str(i) ' ~]=EF_GetPeriod(r1-r0,evs{i}/1000,[wnd(1) wnd(2)],info.Fs,1);']);
    eval(['av' num2str(i) '=EF_XCorr(evs{i},spk,wnd,1,1000);']);
end
xs=wnd(1):1/info.Fs:wnd(2)*2;
s=size(m1);
xs=xs(1:s(2));
% f1=figure;
f2=figure;
for i=1:nn
%     figure(f1);
%     subplot (floor(sqrt(nn))+1,round(sqrt(nn)),i);
%     m=eval(['m' num2str(i)])';
%     plot (xs,m);
%     hold on;
%     plot (xs,mean(m'),'Color','k','LineWidth',2);
%     if ~isempty(patt)
%         hold on;
%         plot(patt(1,:),min(m(:))+(patt(2,:)-min(patt(2,:)))*(max(m(:))-min(m(:))),'r');
%     end
    figure(f2);
    subplot (floor(sqrt(nn))+1,round(sqrt(nn)),i);
    av=eval(['av' num2str(i)]);
    plot (wnd(1)+0.001:.001:wnd(2),MC_smooth(av,1));
    if ~isempty(patt)
        hold on;
        plot(patt(1,:),min(av(:))+(patt(2,:)-min(patt(2,:)))*(max(av(:))-min(av(:))),'r');
    end
end

