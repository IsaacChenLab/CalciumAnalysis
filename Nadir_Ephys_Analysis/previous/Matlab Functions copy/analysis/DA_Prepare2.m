function DA_Prepare2 (fl)

load (fl);
if ~exist('SWEEP_DATA','var') | ~exist('TRI','var')
    disp ('no variables saved');
    return
end
if fl(2)==':'
    nfl=[fl(1:end-4) '_prepared2.mat'];
else
    nfl=['d:\Work\Project 1\Dopamine\' fl(1:end-4) '_prepared2.mat'];
end
s=size(SWEEP_DATA);
for i=2:s(2)
    c1=(SWEEP_DATA(4501,i)-SWEEP_DATA(4000,i))/(TRI(4005,1)-TRI(5000,1));
    SWEEP_DATA(:,i)=SWEEP_DATA(:,i)./c1+TRI(1:5000,2);
end
SWEEP_DATA=SWEEP_DATA(:,2:end);
save (nfl, 'SWEEP_DATA', 'SWEEP_TIMES', 'EVENTS', 'TRI');