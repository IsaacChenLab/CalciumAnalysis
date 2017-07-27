function rez=EF_FiringRate (tms, wnd, res, fq)

% calculates the firing rate based on the spike times provided in tms (ms
% resolution)
% tms - times of spikes
% wnd - window of averaging (sec)
% res - resolution, sliding window
% fq  - frequency used for tms (1000 Hz by default)

if nargin <4
    fq=1000;
end
rez=[];
if res>wnd
    disp ('You selected a resolution larger than the average window. \nThis would result in loss of spikes and inaccurate results.');
    return;
end
n=0;
wnd=wnd*fq;
res=res*fq;
rep=1;

while rep
    rez(end+1,1)=fq*length(find (tms>n & tms<n+wnd))/wnd;
    rez(end,2)=n/fq;
    n=n+res;
    if n>tms(end)
        rep=0;
    end
end