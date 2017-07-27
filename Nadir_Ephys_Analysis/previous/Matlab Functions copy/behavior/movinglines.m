%oriz,for i=0:-0.01:-10, axis([1 5 i i+10]);, pause(0.01);,end,close

ORIZONTAL
n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
t(find(t==0))=10;
t=t-10;
oriz,for i=1:5052, axis([1 5 t(i) t(i)+10]);, pause(0.01);,end,close

ORIZONTAL1
m=2;n=-10:.25:10; 10 epaisseur
n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
t(find(t==0))=10;
t=t-10;
oriz1,for i=1:length(t), axis([1 2 t(i)/2 (t(i)+10)/2]);, pause(0.01);,end,close

ORIZONTAL2
m=2.75;n=-10:.25:10; 10 epaisseur
n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
t(find(t==0))=10;
t=t-10;
oriz2,for i=1:632, axis([1 2 t(i)/2 (t(i)+10)/2]);, pause(0.01);,end,close

VERTICAL
n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
t(find(t==0))=10;
t=t-10;
vert,for i=1:5685, axis([t(i) t(i)+10 0 10]);, pause(0.01);,end,close

ORIZONTAL12
m=1;n=-10:.25:10; 10 epaisseur
n=0;for i=1:0.02:1000, n=n+1;t(n)=i-floor(i/10)*10;,end
t(find(t==0))=10;
t=t-10;
oriz12,for i=1:632, axis([1 2 t(i)/2 (t(i)+10)/2]);, pause(0.01);,end,close

