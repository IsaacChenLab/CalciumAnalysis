function ccc=script_katia
windows=[-500 500];
cats={'c369','c369'}; %,'c370','c370'};
fls={[12 13],[31 32],[21 22],[16 17 18]};
spp=[];
qqq=1;
for k=1:length(cats)
    
    dcat(cats{k});
    sp=cell(1,2);
    groups={[17:24],[25:36]};
    for j=1:2,
        for i=groups{j},
            [s]=MC_loadSpikes(fls{k},i);
            if ~isnan(s)
                sp{j}=[sp{j} s];
            end
        end
    end
    
    for b=1:size(sp{1},2),
        q1=sp{1}(:,b);
        for f=1:size(sp{2},2)
            q2=sp{2}(:,f);
            [c,m,s]=MC_digOnDig(q1,q2,windows);
            ccc(qqq)=sum(c(500:end))-sum(c(1:500));
            qqq=qqq+1;
        end
    end
    [c,m,s]=MC_digOnDig(sum(sp{1},2),sum(sp{2},2),windows);
    subplot(max(length(cats),3),3,(k-1)*3+1);
    xx=windows(1):windows(end);
    cc=c(1:end-1);
    kkk=10;
    cc=sum(reshape(cc,kkk,round(length(cc)/kkk))); 
    xxx=kkk*([1:length(cc)]-round(length(cc)/2));
    cc=cc/mean(cc([1:20]));
    bar(xxx,cc);
    spp=[spp; cc/mean(cc)];
    set(gca,'xlim',[xxx(1) xxx(end)]);
    subplot(max(length(cats),3),3,(k-1)*3+2);
    m=MC_smooth(cc',1);
    plot(xxx,m);
    set(gca,'xlim',[xxx(1) xxx(end)]);
    %set(gca,'ylim',[0.8 1.2]);
    subplot(max(length(cats),3),3,(k-1)*3+3);
    m=MC_smooth(cc',3);
    plot(xxx,m);
    set(gca,'xlim',[xxx(1) xxx(end)]);
    %set(gca,'ylim',[0.8 1.2]);
    clear sp;

end

return;

