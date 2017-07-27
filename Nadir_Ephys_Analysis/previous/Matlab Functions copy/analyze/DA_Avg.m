function [shp se]=DA_Avg (dat,tms,ev,wnd)

s=size(dat);
shp=zeros(s(1),(wnd(2)-wnd(1))*10);
n=0;
for i=1:length(ev)
    p=find(tms>ev(i));
    p=p(1);
    try
        if length(ev)>10000 || s(1)~=1
            shp=shp+dat(:,p+1+wnd(1)*10:p+wnd(2)*10);
        else
            shp(i,:)=dat(:,p+1+wnd(1)*10:p+wnd(2)*10);
            shp(i,:)=shp(i,:)-mean(shp(i,1:0-wnd(1)*10));
        end
        n=n+1;
    catch
    end
end
if length(ev)>10000
    shp=shp/n;
    shp=shp-mean(shp(1:0-wnd(1)*10));
    se=NaN;
elseif s(1)~=1
    se=NaN;
    shp=shp/n;
else
    se=std(shp,0,1)/sqrt(n);
    shp=mean(shp,1);
end