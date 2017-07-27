function MC_plotStim(fls)

groups={[13:16 9:12 5:8 1:4],[21:24 17:20],[25:36],[37:44]};

window=[-250 500];
fl1=fls(1);
[dataStream,totalMs,Sampling_Hz]=MC_openFile(fl1);
[eMat,Times]=MC_loadEvent('Event',fls,1);
eMat=eMat(1:totalMs);

for k=1:length(groups)
    d=1;,
    figure; ca=[]; mn=0; mx=0; d=1;
    for i=groups{k}

        [data,xtime]=MC_getElectrodesInMs(fl1,i);
        [mean_,std_]=MC_anaOnDig(data,eMat,window);
        
        subplot(4,4,d); 
        hold on;
        plot(window(1):window(end),mean_,'k');
        %plot([0 0],[-1000 1000],'-.k');
        set(gca,'xlim',window);
        resp=mean_([1:(abs(window(1))-20) (abs(window(1))+20):end]);
        set(gca,'ylim',[mean(resp)-12*std(resp) mean(resp)+12*std(resp)]);

        d=d+1;

    end
    textglobal(0.5,0.95,num2str(fl1));
    
end

return;

    
    