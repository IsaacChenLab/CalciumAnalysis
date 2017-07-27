function MC_plotSharp1(ind)


if ind==1,
dcat('c369');
fl=[31 32]; fls=[31 32];
groups={[13:16 9:12 5:8 1:4],[21:24 17:20],[25:36]};
ch=15;
elseif ind==2
dcat('c370');
fl=[21 22]; fls=[21 22];
groups={[13:16 9:12 5:8 1:4],[21:24 17:20],[25:36],[37:44]};
ch=14;
elseif ind==3
    dcat('c369');
    fl=[31 32]; fls=[31 32];
    groups={[17:24],[31:34]};
    ch=34;
elseif ind==4
    dcat('c369');
    fl=[12 13]; fls=[12 13];
    groups={[17:24],[33:34]};
    ch=34;
elseif ind==5,
    dcat('c370');
    fl=[21 22]; fls=[21 22];
    groups={[17:24],[33:34]};
    ch=34;
elseif ind==6,
    dcat('c370');
    fl=[21 22]; fls=[21 22];
    groups={[17:24],[33:34]};
    ch=27;
elseif ind==7,
    dcat('c370');
    fl=[16 17 18]; fls=[16 17 18];
    groups={[17:24],[33:34]};
    ch=34;
elseif ind==7,
    dcat('c370');
    fl=[16 17 18]; fls=[16 17 18];
    groups={[17:24],[33:34]};
    ch=27;
end

cdc;
window=[-400 400];
[eMat]=MC_loadEvent('Sharp',fls,ch);
[dataStream,totalMS,Sampling_Hz]=MC_openFile(fls(1));
if length(fl)==1
    xx=find(fls==fl);
    if xx==1
        eMat=eMat(1:totalMS);
    else
        eMat=eMat(totalMS+1 : end);
    end
end
for k=1:length(groups),
    figure; ca=[]; mn=0; mx=0; d=1;
    for i=groups{k}
        data=[];
        for j=1:length(fl)
            [dt,xtim]=MC_getElectrodesInMs(fl(j),i);
            data=[data; dt];
            clear dt;
        end
%         if length(data)>length(eMat)
%             warning('lengths not equal, resampling data');
%             ll1=length(data);
%             ll2=length(eMat);
%             data=decimate(data,round(ll1/ll2));
%         end
        [mean_,std_]=MC_anaOnDig(data,eMat,window);
        clear data;
        ca{d}=subplot(4,4,d); 
        hold on;
        plot(window(1):window(end),mean_,'k');
        plot([0 0],[-1000 1000],'-.k');
        set(gca,'xlim',window);
        mn=min(mn,min(mean_));
        mx=max(mx,max(mean_));
        for ii=1:d,
            %axes(ca{d}); hold on;
            set(ca{ii},'ylim',[mn mx]);
            hold off;
        end
        d=d+1; 
    end
    drawnow;
end

return;