function rez= XCorr (ev_ref, ev2, wnd, bin)

if nargin<4
    error ('not enough arguments');
end
frez=[];
for i=1:length(ev_ref)
    tev2=ev2-ev_ref(i);
    frez=[frez tev2(tev2>wnd(1) & tev2<wnd(2))];
end
xs=wnd(1):bin:wnd(2);
[rez b]=hist(frez,xs);
rez=rez./length(ev_ref);
end