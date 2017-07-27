function []=DA_PlotEvents (EVENTS,SWEEP_DATA,SWEEP_TIMES,sh)

s=size(SWEEP_DATA);
for i=1:s(2)
    dat(i)=dot(SWEEP_DATA(:,i),sh);
end
clear SWEEP_DATA;

s=size(EVENTS);
ev=sparse(zeros(1,s(2)));
ev(10000:100:s(2))=1;
[ref ~]=DA_Avg (dat,SWEEP_TIMES,ev,[-5 10]);

nms={'LKS'; 'RWD'; 'CSP'; 'CSM'};
col={'k'; 'g'; 'b'; 'r'};

rez=zeros(4,150);
err=rez;
f1=figure;
hold on
f2=figure;
hold on

for i=1:4
    [rez(i,:) err(i,:)]=DA_Avg (dat,SWEEP_TIMES,EVENTS(i,:),[-5 10]);
    figure(f1);
    plot (-4.9:0.1:10,rez(i,:),[col{i} '-']);
    plot (-4.9:0.1:10,rez(i,:)+err(i,:),[col{i} '--']);
    plot (-4.9:0.1:10,rez(i,:)-err(i,:),[col{i} '--']);
    figure(f2);
    plot (-4.9:0.1:10,rez(i,:)-ref,[col{i} '-']);
    plot (-4.9:0.1:10,rez(i,:)+err(i,:)-ref,[col{i} '--']);
    plot (-4.9:0.1:10,rez(i,:)-err(i,:)-ref,[col{i} '--']);
end


