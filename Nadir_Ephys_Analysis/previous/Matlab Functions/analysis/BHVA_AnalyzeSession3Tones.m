function rez=BHVA_AnalyzeSession3Tones (fn)

try
    load (fn);
catch
    error ('file could not be loaded');
end
if ~(exist('events','var') && exist ('licks','var') && exist ('i','var'))
    error ('the loaded session does not contain appropriate variables');
end
try
EVENTS{1}=events{1,1};
catch
end
try
EVENTS{2}=events{1,2};
catch
end
try
EVENTS{3}=events{1,3};
catch
end
try
EVENTS{4}=events{1,4};
catch
end
try
LICKS=licks(1,:);
catch
end

rt=Rate(licks,10)/10;
% EVENTS=events;
f=figure;
plot (rt);
hold on;
[x1 ~]=ginput(1);
line([x1 x1],[min(rt) max(rt)]);
[x2 ~]=ginput(1);
close (f);
MX=round(10*max(x1,x2));
MN=round(10*min(x1,x2));

for ii=1:4
    try
        ev=eval(['EVENTS{' num2str(ii) '}']);
        EV=ev(ev<MX);
        EV=EV(EV>MN);
        eval (['EVENTS{' num2str(ii) '}=EV;']);
    catch
    end
end

for ii=1:4
    try
        eval (['ev' num2str(ii) '=MC_smooth(XCorr (EVENTS{' num2str(ii) '},LICKS,[-5 10],.1)'',0.5);']);
        eval(['rez.evp' num2str(ii) '=ev' num2str(ii) ';']);
    catch
    end
end

mx=0;
for ii=1:4
    try
        mx=max(mx,max(eval(['ev' num2str(ii)])));
    catch
    end
end
figure

legend={'REWARD'; 'CS1'; 'CS2';'CS3'};
for ii=1:4
    try
        subplot (2,2,ii);
        plot (-5:.1:10,eval(['ev' num2str(ii)]));
        axis ([-5 10 0 mx]);
        title (legend{ii});
    catch
    end
end



% plot (0:10:LICKS(end),rt);
% hold on
% line([x1*10 x1*10],[min(rt) max(rt)],'Color','r');
% line([x2*10 x2*10],[min(rt) max(rt)],'Color','r');
% axis ([0 LICKS(end) 0 max(rt)]);
% title ('Licking rate');
% 
%%%% percent hits/miss
ev=[];
for ii=1:length(EVENTS{1})
    if ~isempty(find(LICKS>EVENTS{1}(ii) & LICKS<EVENTS{1}(ii)+1))
        ev(end+1)=EVENTS{1}(ii);
    end
end
EVENTS{1}=ev;
for ii=2:4
    try
        ev=eval(['EVENTS{' num2str(ii) '}']);
        n=0;
        for j=1:length(ev)
            if ~isempty (find(EVENTS{1}>ev(j) & EVENTS{1}<ev(j)+5))
                n=n+1;
            end
        end
        eval(['rez.ev' num2str(ii) '=n/length(ev);']);
        eval (['n' num2str(ii-1) '=n;']);
    catch
    end
end
rez.cor=(n1+n2)/(length(EVENTS{2})+length(EVENTS{3}));
rez.l2=sum(ev2(51:62));
rez.l3=sum(ev3(51:62));
rez.l4=sum(ev4(51:62));
lp=sum(ev2(51:62));
blp=sum(ev2(41:50));
lm=sum(ev3(51:62));
blm=sum(ev3(41:50));
rez.dis=(lp-lm)/(lp+lm);
rez.percent_lp=lp/blp;
rez.percent_lm=lm/blm;
