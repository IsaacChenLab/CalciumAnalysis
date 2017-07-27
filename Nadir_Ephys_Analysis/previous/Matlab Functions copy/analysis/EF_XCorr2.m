function rez=EF_XCorr2 (ref, dat, wnd, bin, fs)

wnd=wnd*fs;
rez=zeros(wnd(2)-wnd(1),length(ref));
for i=1:length(ref)
    d=dat(dat>ref(i)+wnd(1))-ref(i)-wnd(1)+1;
    if ~isempty(d(d<wnd(2)-wnd(1)))
        rez(round(d(d<wnd(2)-wnd(1))),i)=1;
    end
end
rez=rez/length(ref);
if bin>1
    rez=MC_bin(rez',bin)';
end