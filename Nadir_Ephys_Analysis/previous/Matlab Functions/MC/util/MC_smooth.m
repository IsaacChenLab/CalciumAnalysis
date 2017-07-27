function [vout,win]=MC_smooth(vinp,width)

if width==0
    vout=vinp;
    wind=[0];
    return;
end

%win=normpdf(-2*width:width:2*width,0,width); 
%win=win/sum(win);

l=min(round(size(vinp,1)/2),ceil(2*width));
win=normpdf(-l:l,0,width); 
win=win/sum(win);

L1=size(vinp,1);
L2=ceil(width*3);
tvinp=[vinp(L2:-1:1,:); vinp; vinp(L1:-1:(L1-L2),:)];
tvout=filtfilt(win,1,tvinp);
vout=tvout(L2:(L2+L1-1),:);

return;

