scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);
drawnow;
m=0; %for horiz, alfa=0; m=tg(alfa=0)=sin0/cos0=0/1=0; if alfa=90 (vert), m=tg(alfa=90)=sin90/cos90=1/0=infinit  
hold on;
for n=-10:0.5:10
for i=1:5
    x1(i)=i;
    y1(i)=m*x1(i)+n;
end
plot(y1,'w','LineWidth',10);subplot('position',[0 0 1 1]);hold on;
end
%axis square;
axis fill;
axis off;