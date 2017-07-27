function xpdf=MC_histToX(mat,norm)

nx = size(mat,1); 
ny = size(mat,2);  
win=[-round(nx/2) round(nx/2)];
x_axis = linspace(win(1),win(2),nx);
y_axis = linspace(win(1),win(2),ny);

[X,Y] = meshgrid(x_axis, y_axis);

X=X(:);
Y=Y(:);
m=mat(:);
if nargin==2 & norm
    m=m-min(m);
end

xpdf=nan(nansum(m),2);
k=1;
for i=1:length(X),
    xpdf(k:k+m(i)-1,1)=X(i)*ones(m(i),1)+rand(m(i),1)-0.5;
    xpdf(k:k+m(i)-1,2)=Y(i)*ones(m(i),1)+rand(m(i),1)-0.5;
    k=k+m(i);
end

return;


