% s=size(dat);
% cf=2;
% m=mean(dat,2);
% for i=1:s(2)
%     vl=dat(:,i).*COEFF(:,cf);
% %     plot (TRI,vl);
%     plot (TRI,dat(:,i)-m);
%     axis ([min(TRI) max(TRI) -1 1]);
%     pause;
% end
% 
% 
% %%%%%%%%%%%%%%%%
% for j=1:150
% s=zeros(2752,1);
% n=0;
% for i=3+j:150:4350
%     n=n+1;
%     s=s+NcAcc_SWEEP_DATA(1:2752,i);
% end
% s=s./n;
% plot (s)
% pause;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%
% for i=1:10
%     s=zeros(150,1);
%     n=0;
%     for j=1:150:3290-150
%         n=n+1;
%         s=s+aSCORE(j:j+149,i);
%     end
%     s=s./n;
%     plot (s);
%     pause
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% s=size(dat);
% cf=4;
% m=mean(dat,2);
% for i=1:s(2)
%     vl=dat(:,i);
% %     plot (TRI,vl);
%     plot (TRI(1:len,1),dat(:,i));
%     axis ([0.2 1 2 10]);
%     pause(0.1);
% end


%%%%%%%%%%%%%% average by xxs sec
xxs=15;
wh=800;
m=zeros(xxs*20,1);
f=1;
pos=1;
n=0;
while f
    val=SWEEP_TIMES(pos)+xxs;
    ff=find (SWEEP_TIMES<val);
    ppos=ff(end)+1;
    if ppos>length(SWEEP_TIMES)
        f=0;
    else
        ex=pos:ppos;
        m(1:length(ex))=m(1:length(ex))+dat(wh,ex)';
        n=n+1;
    end
    pos=ppos+1;
end
m=m./n;
plot (m(1:xxs*10)./mean(m(1:xxs*10)),'r')