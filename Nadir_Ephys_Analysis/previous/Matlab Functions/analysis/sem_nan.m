function rez=sem_nan(mat)

s=size(mat);
rez=zeros(s(2),1);

for i=1:s(2)
    rez(i)=std(mat(~isnan(mat(:,i)),i))/sqrt(length(mat(~isnan(mat(:,i)),i)));
end