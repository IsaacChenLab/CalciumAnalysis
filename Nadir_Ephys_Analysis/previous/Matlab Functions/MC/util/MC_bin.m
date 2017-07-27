function d1=MC_bin(data,bin)

if bin==1,
    d1=single(data); 
    clear data;
    return;
end

% if size(data,2)>size(data,1)
%     data=data';
%     rev=1;
% else
%     rev=0;
% end

if size(data,2)>1
    rev=1;
    data=data';
else
    rev=0;
end

if mod(size(data,1),bin)>0
    data=data(1:bin*floor(size(data,1)/bin),:);
    warning('data truncated because of bin size');
end
[N,M]=size(data);
CHUNKS=10000;
if mod(CHUNKS,bin)>0
    while mod(CHUNKS,bin)>0
        CHUNKS=CHUNKS+1;
    end
    warning(['CHUNKS increased to ' num2str(CHUNKS)]);
end

periods=floor(M/CHUNKS); last=mod(M,CHUNKS);
d1=zeros(N/bin,M,'single');
for i=1:periods,
    d2=reshape(data(:,[1:CHUNKS]+(i-1)*CHUNKS),bin,CHUNKS*N/bin);
    d2=single(sum(d2));
    d1(:,[1:CHUNKS]+(i-1)*CHUNKS)=reshape(d2,N/bin,CHUNKS);
end
d2=reshape(data(:,[periods*CHUNKS+1:periods*CHUNKS+last]),bin,last*N/bin);
d2=single(sum(d2));
d1(:,[periods*CHUNKS+1:periods*CHUNKS+last])=reshape(d2,N/bin,last);
clear data d2;

d1=double(d1);

if rev
    d1=d1';
end

return;
