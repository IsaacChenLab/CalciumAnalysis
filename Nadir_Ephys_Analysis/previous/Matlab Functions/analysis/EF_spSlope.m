function inxs=EF_spSlope(wv)

av=mean(wv,2);
f=figure;
rep=1;
while rep
    hold off;
    plot (av,'k');
    hold on;
    x=zeros(1,3);
    for i=1:3
        [x(i), ~]=ginput(1);
        x(i)=min(73,max(1,floor(x(i))));
        plot (x(i),av(x(i)),'r.','MarkerSize',30);
    end
    q=questdlg ('Accept these values?', 'Points set.','Yes','No','No');
    if strcmp(q,'Yes')
        rep=0;
    end
end
close (f);
x=sort(x);
for i=1:size(wv,2)
    sl(1,i)=wv(x(1),i)-wv(x(2),i);
    sl(2,i)=wv(x(2),i)-wv(x(3),i);
end
sl=abs(sl);


val1= sl(1,:);
val2= sl(2,:);
ind2=1;

f1=figure; hold on;

while ind2,
    plot(val1,val2,'.k','markersize',5);
%     set(gca,'xtick',[]); set(gca,'ytick',[]);
    axis square;

    ind1=1;
    d=1;
    points=[];
    while ind1,
        [x,y]=ginput(1);
        points(d,:)=[x y];
        if d>1
            h=line([points(d-1,1) points(d,1)],[points(d-1,2) points(d,2)]);
            set(h,'color','b');
        end
        if d>3
            xxx=MC_spFindCross(points(1:d-1,1),points(1:d-1,2),points(d-1:d,1),points(d-1:d,2));
            if ~isempty(xxx)
                ind1=0;
                break;
            end
        end
        d=d+1;
    end

    IN=inpolygon(val1,val2,points(:,1),points(:,2));
    plot(val1(IN),val2(IN),'.b','markersize',5);
    s=questdlg('Accept?','','Yes','No','Yes');
    if strcmp(s,'Yes'),
        break;
        ind2=0;
    elseif strcmp(s,'No')
        cla;
    elseif isempty(s)
        error('');
    end
    
end
close(f1);
inxs{1}=find(IN);
inxs{2}=find(~IN);