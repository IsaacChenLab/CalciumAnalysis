function [cc,cc_smooth,lags,between]=MC_rawCC(spikeMat,cc_width)

if nargin<2
    cc_width=500; % 500 ms
end

lags=[-cc_width:cc_width];

[nSamples,nCells]=size(spikeMat);

cc=[];
d1=1;
for i=1:nCells,
    s1=spikeMat(:,i);
    xs=find(s1);
    for j=i:nCells,
        between(d1,:)=[i j];
        cc=[cc zeros(2*cc_width+1,1)];
        s2=full(spikeMat(:,j));
        d2=1;
        for k=1:length(xs),
           q=[xs(k)-cc_width : xs(k)+cc_width];
           if q(1)>0 & q(end)<length(s2)
               cc(:,d1)=cc(:,d1)+s2(q);
               d2=d2+1;
           end
        end
        if i==j, % auto
            cc(cc_width+1,d1)=0;
        end
        cc(:,d1)=cc(:,d1)/d2;
        d1=d1+1;
    end
end

cc_smooth=MC_smooth(cc,round(cc_width/40));

return;

        