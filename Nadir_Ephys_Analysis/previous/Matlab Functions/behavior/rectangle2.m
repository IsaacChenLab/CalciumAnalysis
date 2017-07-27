function H=rectangle(a,b,style)
X=[-a a a -a];
Y=[b b -b -b];
fill(X,Y,style); hold on;
H=plot(X,Y,style);
axis equal;
axis off;
