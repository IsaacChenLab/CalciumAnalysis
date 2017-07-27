function H=circle(center,radius,NOP,style)
if (nargin <3),
 error('Please see help for INPUT DATA.');
elseif (nargin==3)
    style='b-';
end;
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);fill(X,Y,style); hold on;
H=plot(X,Y,style);
axis equal;
axis square;
axis off;
