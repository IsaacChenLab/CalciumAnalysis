function [counts,me,mean_,std_]=MC_digOnDig(event1,event2,wind)

%event1=full(event1);
%event2=full(event2);
if length(event1)~=length(event2)
    warning('lengths not equal, truncating according to the shorter one');
    ll=min(length(event1),length(event2));
    event1=event1(1:ll);
    event2=event2(1:ll);
end
if nargin<3
    wind=-500:500;
elseif length(wind)==2
    wind=wind(1):wind(2);
elseif length(wind)==1
    wind=-wind:wind;
end

jx=find(event2);
if ~isempty(jx)
    clear event2;
    jx=jx(find( (jx+wind(1))>0 & (jx+wind(end))<=length(event1) ));
    jx=single(jx); wind=single(wind);
    CHUNKS=10000;
    periods=floor(length(jx)/CHUNKS); last=mod(length(jx),CHUNKS);
    me=zeros(length(jx),length(wind),'uint8');
    for i=1:periods,
        mx=repmat(jx([1:CHUNKS]+(i-1)*CHUNKS),1,length(wind)) + repmat(wind,CHUNKS,1);
        me([1:CHUNKS]+(i-1)*CHUNKS,:)=uint8(full(event1(mx)));
    end
    mx=repmat(jx([periods*CHUNKS+1:periods*CHUNKS+last]),1,length(wind)) + repmat(wind,last,1);
    me(periods*CHUNKS+1:periods*CHUNKS+last,:)=uint8(full(event1(mx)));
    clear jx mx event1;
else
    me=nan(1,length(wind));
    warning('no ref spikes');
end

me=double(me);
counts=(sum(me,1));
if nargout>2
    mean_=(mean(me,1));
    std_=std((me),0,1);
end


return;






