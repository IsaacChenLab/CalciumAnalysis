function [mout,p]=MC_smooth(minp,width)

[x,y]=meshgrid(-2*width:2*width);
for i=1:size(x,1)
    for j=1:size(x,2),
        p(i,j)=mvgauss([x(i,j) y(i,j)],[0 0],[width 0; 0 width]);
    end
end

p=p/sum(p(:));

mout=filter2(p,minp);

return;

