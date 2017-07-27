function [vv,inx,shifts]=MC_spRealign(v,Fs)
% m*n, samples * incidents 

if nargin<2
    Fs=25e3;
end

v=v';

Minus=round(Fs/1e3);
Plus=round(2*Fs/1e3);

if isempty(v)
    vv=[];
    return;
end

v=v';
[mn,xmn]=min(v);


d=1; inx=[];
for i=1:size(v,2)
    if xmn(i)-Minus<= 0 | xmn(i)+Plus > size(v,1),
        ;
    else
        vv(:,i)=v(xmn(i)-Minus : xmn(i)+Plus-1,i);
        inx(d)=i;
        shifts(d)=xmn(i)-Plus;
        d=d+1;
    end
end

if isempty(inx)
    vv=[];
    return;
end

vv=vv(:,inx);

return;


