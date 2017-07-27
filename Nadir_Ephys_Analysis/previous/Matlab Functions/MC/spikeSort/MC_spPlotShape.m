function MC_spPlotShape(cwv,cf,mnn,mxx,c)

if isempty(cwv),
    return;
end

if nargin<3
    mnn=min(cwv(:));
    mxx=max(cwv(:));
end
if nargin<5
    c='k';
end
if nargin>1
    axes(cf);
end
hold on;
plot(cwv,'color',c);
mn=mean(cwv'); 
st=std(cwv');
nSamples=size(cwv,1);
set(gca,'XLim',[1 nSamples]);
set(gca,'XTick',[1 round(nSamples/3) round(2*nSamples/3) nSamples]);
set(gca,'XTickLabel',{'-1 ms','0','+1ms','+2ms'});
set(gca,'YLim',[min(cwv(:))-range(cwv(:))/10 max(cwv(:))+range(cwv(:))/10]);
% set(gca,'YLim',[min(cwv(:)) max(cwv(:))]);
ylabel('\mu volt');
plot(mn,'color',[0.5 0.5 0.5],'LineWidth',3);
[ttt,x1]=min(mn); [ttt,x2]=max(mn(26:end)); x2=x2+25; x=[12 x1 x2];
errorbar(x,mn(x),st(x),'.','LineWidth',3,'color',[0.5 0.5 0.5]);
% set(gca,'ylim',[mnn mxx]);
text(nSamples-20,mxx-round((mxx-mnn)/20),['n=' num2str(size(cwv,2))]);

return;


