function H=triangle(a,style)
X=[0 -a/2 a/2];
Y=[0.865*a 0 0];
fill(X,Y,style); hold on;
H=plot(X,Y,style);
axis off;
