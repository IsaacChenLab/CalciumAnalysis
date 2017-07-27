function [Inxs]=MC_spAmplitude(waveForms)

if isempty(waveForms)
    error('no spikes provided');
end

mins=min(waveForms(:,:));
%maxs=max(waveForms(:,:));

f1=figure;
xaxis=1:length(mins);
%subplot(2,1,1);
plot(xaxis,mins,'k.'); set(gca,'xlim',[1 length(mins)]);
%subplot(2,1,2);
%plot(xaxis,maxs,'k.');

c={'b','r','g','y'};
d=1;
ind=1;
while ind,
    figure(f1); hold on;
    y=[0 0];
    [x,y]=ginput(2);
    if ~isempty(x),
        y=sort(y);
        if d>1 & y(2)>ysave
            y(2)=ysave;
        end
        xx=find(mins>=y(1) & mins<y(2));
        plot(xaxis(xx),mins(xx),[c{d} '.']);
        ysave=y(1);
        Inxs{d}=xx;
        d=d+1; 
    else
        ind=0;
    end
end

close(f1);

return;


