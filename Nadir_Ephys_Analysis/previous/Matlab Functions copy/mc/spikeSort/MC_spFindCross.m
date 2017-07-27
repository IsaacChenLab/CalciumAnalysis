function xxx=MC_spFindCross(Xs,Ys,x,y)

[nx,ny]=size(Xs);
Xs=Xs(:);
Ys=Ys(:);
Xs2=circshift(Xs,[-1 1]);
Ys2=circshift(Ys,[-1 1]);
% find line formula for all pairs of points in the shapes
% y=Ax+B;
As=(Ys2-Ys) ./ (Xs2-Xs);
Bs=Ys-(As.*Xs);
% y=ax+b for the line [x y]
a=(y(2)-y(1))/(x(2)-x(1));
b=y(1)-(a*x(1));
% find x intersect of line and all lines from shapes, by ax+b=Ax+B
Xintersect=(Bs-b)./(a-As);

x=sort(x);
xxx=find(Xintersect>=x(1) & Xintersect<=x(2) & Xintersect>Xs-eps & Xintersect<Xs2+eps);
md=mod(xxx,(nx));
xxx=xxx(md>0); % avoid edges 
xxx=ceil(xxx/(nx));
xxx=unique(xxx);

return;
