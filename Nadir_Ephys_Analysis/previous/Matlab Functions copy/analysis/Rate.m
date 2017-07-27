function rez=Rate (ev, bin)

xs=0:bin:max(ev);
[rez b]=hist(ev,xs);