function rez=BHV_AnalyzeTone3 (fl, freqq)

load (fl)
if ~strcmp(Proc,'T3')
    error ('Incorrect session type');
end

cspt=events{2};
csmt=events{3};
freq=events{11};
cs=[cspt csmt];
[~, idx]=sort (cs);
csmf=freq(idx>length(cspt));
if strcmp(questdlg('Laser used?','Question','Yes','No','No'),'No')
    laser=0;
else
    laser=1;
end

e1=XCorr (cspt, licks, [-3 7], .1);
e2=XCorr (csmt, licks, [-3 7], .1);
e3=XCorr (cspt, events{1}, [-3 7], .1);
f=figure;
hold on;
plot (-3:.1:7,e1,'b');
plot (-3:.1:7,e2,'k');
plot (-3:.1:7,e3,'r');
[x, ~, ~]=ginput(2);
close (f)
x=[min(x) max(x)];
tmp=find(freqq==800);
% rez.freq=800-tmp:800+tmp;
rez.freq=freqq;
rez.licks=zeros(1,length(rez.freq));
nn=rez.licks;
for i=1:length(csmt)
p=round(csmf(i)-800+tmp);
nn(p)=nn(p)+1;
rez.licks(p)=rez.licks(p)+length(find(licks>csmt(i)+x(1) & licks<csmt(i)+x(2)));
end
figure;
plot (rez.freq,MC_smooth(rez.licks',10));