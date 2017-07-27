function EF_FixEvents(id,rt)

if nargin<2
    rt='H:';
end

if ~exist([rt '\To Matlab\' id '\digi01.mat'],'file') || ~exist([rt '\BHV\' id '.mat'],'file')
    warndlg('No such files.');
    return
end

load([rt '\To Matlab\' id '\digi01.mat']);
load ([rt '\BHV\' id '.mat']);
tms=0:1/info.Fs:length(digi)/info.Fs;
df=diff(digi);
f=find(df<0);
f=tms(f);
e=events{1};
ref=1;
rep=1;
gf=figure;
while rep
    figure(gf);
    hold off;
    sft=e(1)-f(ref);
    plot (f+sft,0,'rx');
    hold on
    plot (e,0,'bo');
    [p pp but]=ginput(1);
    if isempty(p)
        rep=0;
    elseif but==1 
        ref=ref-1;
    elseif but==3 && ref<length(f)
        ref=ref+1;
    end
end
close (gf);
sft=f(ref)-e(1);
s=size(events);
EVENTS=[];
EVENTS{1}=licks+sft;
EVENTS{1}=EVENTS{1}(EVENTS{1}>0);
for i=1:s(2)
    if ~isempty(events{i})
        EVENTS{i+1}=events{i}+sft;
        EVENTS{i+1}=EVENTS{i+1}(EVENTS{i+1}>0);
    end
end
save ([rt '\To Matlab\' id '\DIG.mat'],'EVENTS');