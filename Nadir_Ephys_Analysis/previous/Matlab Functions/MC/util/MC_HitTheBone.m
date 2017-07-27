function MC_HitTheBone

[cats,dates,files]=MC_getDataForSpike;
q=1;
for i=1:length(cats),
    dcat(cats{i},dates{i});
    if strcmp(cats{i},'c375') | strcmp(cats{i},'c378'),
        d(q)=str2num(dates{i});
        for j=1:16
            [SpikeMat,SpikeGrade]=MC_loadSpikes(files{i},j);
            if isempty(SpikeMat) | isnan(SpikeMat)
                n(j,q)=0;
            else
                n(j,q)=length(SpikeGrade);
            end
        end
        q=q+1;
    end
end

dc=repmat(3*[0:size(n,1)-1]',1,length(d));
n=n+dc;
figure;
x1=find(d<700);
x2=find(d>700);

subplot(1,2,1); hold on;
plot(d(x1),n(:,x1));
plot(d(x1),dc(:,x1),'k--');

subplot(1,2,2); hold on;
plot(d(x2),n(:,x2));
plot(d(x2),dc(:,x2),'k--');

return;