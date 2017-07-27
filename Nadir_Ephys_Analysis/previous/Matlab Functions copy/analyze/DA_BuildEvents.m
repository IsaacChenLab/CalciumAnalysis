function DA_BuildEvents (fl)

load (fl,'EVENTS');
s=size(EVENTS);
answer=inputdlg(['There are ' num2str(s(1)) ' events recorded. Monitor:'], ' ',1,{'1'});
try
    a=str2num(answer);
catch
    return;
end
if a<1 || a>s(1)
    return;
end
new_ev=EVENTS(a,:);

evs=EVENTS;
evs(a,:)=0;
flt=DA_MakeFLT (evs,new_ev);
EVENTS(end+1,:)=EVENTS(a,:).*flt;
save (fl,'EVENTS','-append');