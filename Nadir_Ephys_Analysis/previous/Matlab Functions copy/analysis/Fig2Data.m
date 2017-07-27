function rez=Fig2Data (fig)

figure(fig);
h = gcf;
axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children');
objTypes = get(dataObjs, 'Type');
rez.x = get(dataObjs, 'XData'); 
rez.y = get(dataObjs, 'YData');
rez.z = get(dataObjs, 'ZData');
rez.type=objTypes;