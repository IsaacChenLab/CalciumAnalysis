function MC_showDataAll(file,nch,prind,notch60)
G=8;
Times={[40 46]};%,[35 36],[35.5 35.6],[35.55 35.575]};
if nargin<2
    nch=1:MC_NumElecInFile(file);
end
if nargin<3
    prind=0;
end
if nargin<4
    notch60=0;
end

[dataStream,totalSeconds,Hz]=MC_openFile(file);
periods=ceil(length(nch)/G);
last=mod(length(nch),G);

for i=1:periods,
    if last>0 & i==periods
        ch=nch([0+[1:last]+(i-1)*G]);
    else
        ch=nch([0+[1:G]+(i-1)*G]);
    end
    [a1,t1]=MC_getElectrodes(file,ch,Times{1});
    if notch60==1
        a1=MC_notch60(a1,Hz);
    end
    for j=1:length(Times),
        f1=figure;
        xx=find(t1>=Times{j}(1) & t1<=Times{j}(2));
        for k=1:length(ch)
            subplot(4,2,k);
            d=a1(xx,k);
            plot(t1(xx),d,'k');
            title(['channel ' num2str(ch(k))]);
            set(gca,'xlim',[t1(xx(1)) t1(xx(end))]);
            set(gca,'xtick',[t1(xx(1)) t1(xx(1))+(t1(xx(end))-t1(xx(1)))/2 t1(xx(end))]);
            if ~isnan(range(d)) & range(d)>0
                set(gca,'ylim',[min(d)-range(d)/5 max(d)+range(d)/5+eps]);
            end
        end
        textglobal(0.5,0.95,sprintf('%2.3f - %2.3f sec',Times{j}(1),Times{j}(end)));
        if prind
            cdc;
            cd('Figures');
            print(f1,'-depsc2',[num2str(file) '_' num2str(ch(1)) '_' num2str(ch(end)) '_' num2str(j)]);
        end
    end
end 

return;

        

