function rez = EF_Convert2ZScore (data, dim)

if dim==1
    data=data';
end

m=mean(data,1);
s=std(data,1);
rez=data;
rez(:)=0;
for i=1:size(data,2)
    rez(:,i)=(data(:,i)-m(i))/s(i);
end