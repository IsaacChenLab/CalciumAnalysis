function [mean_,std_,me]=MC_anaOnDig(data,event,window)

event=full(event);
if length(data)>length(event)
    warning('lengths not equal, resampling data');
    ll1=length(data);
    ll2=length(event);
    data=resample(data,ll2,ll1);
end
if nargin<3
    window=-500:500;
elseif length(window)==2
    window=window(1):window(2);
elseif length(window)==1
    window=-window:window;
end

jx=find(event);
jx=jx(find( (jx+window(1))>0 & (jx+window(end))<=length(data) ));
mx=repmat(jx,1,length(window)) + repmat(window,length(jx),1);
me=data(mx);

mean_=mean(me);
std_=std(me);

return;
