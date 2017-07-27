function CPP_Analyse (ss)

try
    load (ss)
    expl=[White.ExplorationCounts Gray.ExplorationCounts Black.ExplorationCounts]';
    entr=[White.EntranceCounts Gray.EntranceCounts Black.EntranceCounts]';
    tms=[White.ZoneTime Gray.ZoneTime Black.ZoneTime]';
    activ=[White.ActivityCounts Gray.ActivityCounts Black.ActivityCounts]';
    mvm=[White.MovementCounts Gray.MovementCounts Black.MovementCounts]';
catch
    error ('File or variables missing');
end

s=size(expl);
for i=3:s(2)
    expl(:,i)=expl(:,i)+expl(:,i-1);
    entr(:,i)=entr(:,i)+entr(:,i-1);
    tms(:,i)=tms(:,i)+tms(:,i-1);
    activ(:,i)=activ(:,i)+activ(:,i-1);
    mvm(:,i)=mvm(:,i)+mvm(:,i-1);
end
figure ('Name',ss)
subplot (3,2,1);
hold on
plot (expl(1,2:end),'r');
plot (expl(2,2:end),'y');
plot (expl(3,2:end),'k');
title ('ExplorationCounts')
try
axis ([1 s(2)-1 0 max(expl(:,1))]);
catch end

subplot (3,2,2);
hold on
plot (entr(1,2:end),'r');
plot (entr(2,2:end),'y');
plot (entr(3,2:end),'k');
title ('EntranceCounts')
try
axis ([1 s(2)-1 0 max(entr(:,1))]);
catch end


subplot (3,2,3);
hold on
plot (tms(1,2:end),'r');
plot (tms(2,2:end),'y');
plot (tms(3,2:end),'k');
title ('ZoneTime')
try
axis ([1 s(2)-1 0 max(tms(:,1))]);
catch end


subplot (3,2,4);
hold on
plot (activ(1,2:end),'r');
plot (activ(2,2:end),'y');
plot (activ(3,2:end),'k');
title ('ActivityCounts')
try
axis ([1 s(2)-1 0 max(activ(:,1))]);
catch end

subplot (3,2,5);
hold on
plot (mvm(1,2:end),'r');
plot (mvm(2,2:end),'y');
plot (mvm(3,2:end),'k');
title ('MovementCounts')
try
axis ([1 s(2)-1 0 max(mvm(:,1))]);
catch end


subplot (3,2,6);
hold on
plot (mvm(1,2:end)./tms(1,2:end),'r');
plot (mvm(2,2:end)./tms(2,2:end),'y');
plot (mvm(3,2:end)./tms(3,2:end),'k');
title ('Normalized Movement')
