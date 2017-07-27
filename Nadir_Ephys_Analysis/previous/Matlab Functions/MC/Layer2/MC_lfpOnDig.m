function [mean_,std_,me]=MC_anaOnDig(data,event,window)

event=full(event);
if nargin<3
    window=-500/2:500/2;
elseif length(window)==2
    window=round(window(1)/2):round(window(2)/2);
elseif length(window)==1
    window=-round(window/2):round(window/2);
end

jx=find(event); 
jx=round(jx/2);
jx=jx(find( (jx+window(1))>0 & (jx+window(end))<=length(data) ));
mx=repmat(jx,1,length(window)) + repmat(window,length(jx),1);
me=data(mx);

mean_=mean(me);
std_=std(me);

return;
