function rez=EF_SPK2LFP (spk,leng,jit)

rez=zeros(1,round(leng*1000));
v=pdf('norm',[-4*jit:4*jit],0,jit);
for i=1:length(spk)
    p=round(spk(i));
    rez(p:p+length(v)-1)=rez(p:p+length(v)-1)+v;
end
rez=rez(ceil(length(v)/2):end);