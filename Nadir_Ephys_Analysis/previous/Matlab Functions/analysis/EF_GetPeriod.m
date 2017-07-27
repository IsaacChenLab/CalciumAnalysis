function [mat tms]=EF_GetPeriod (data, events, wnd, Fs, sr)

% data - vector of continuous recording to be segmented
% events - time points (in sec) used to center the segments
% wnd - [pre post] time windows (in sec) to retain pre-event and post-event
% Fs - sampling frequency of data
% sr - 1 means select all events, 0 means choose manually

if nargin<5
    sr=[];
end
evs=round(events*Fs);
win=round(wnd*Fs);
win=[min(win) max(win)];
mat=zeros(length(evs),win(2)-win(1)+1);
nw=[];
for i=1:length(evs)
    try
        mat(i,:)=data(evs(i)+win(1):evs(i)+win(2));
    catch
        nw(end+1)=i;
    end
end
mat=mat(setdiff(1:length(evs),nw),:);

if isempty (sr)
    q=questdlg ('Do you want to select the events?','Question','Yes','No','No');
elseif sr
    q='No';
else
    q='Yes';
end
    
if strcmp(q,'No')
    tms=events;
    return
end
n=[];
f=figure;
hold off;
for i=1:length(events)
    plot (mat(i,:));
    [x xx b]=ginput(1);
    if b==3
        n(end+1)=i;
    end
end
mat=mat(setdiff(1:length(events),n),:);
tms=events(setdiff(1:length(events),n));
end
