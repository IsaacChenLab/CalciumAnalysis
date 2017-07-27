scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);
drawnow;
m=2.75; 
hold on;
for n=-10:.25:10
for i=1:5
    x1(i)=i;
    y1(i)=m*x1(i)+n;
end
plot(y1,'w','LineWidth',10);subplot('position',[0 0 1 1]);hold on;
end
%axis square;
axis fill;
axis off;