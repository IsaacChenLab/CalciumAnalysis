scrsz = get(0,'ScreenSize');
f1=figure('Position',scrsz,'color',[0 0 0]);
drawnow;
for a=-10:.3125:10
    x=[a a a a a a a a a a a];
    y=[0 1 2 3 4 5 6 7 8 9 10];
    hold on;
    plot(x,y,'w','LineWidth',10);subplot('position',[0 0 1 1]);hold on;
end
%axis square;
%axis equal;
axis fill;
axis off;