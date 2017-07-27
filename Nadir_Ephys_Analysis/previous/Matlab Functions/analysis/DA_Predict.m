function [rez, rezc, nm, nmc]=DA_Predict (fl, pnts, mod)

% try to predict whether the animal will lick for reward or not based on
% the dopamine signal
%
%fl - prepared file to load
%pnts - data points following CS to be used for prediction
%mod - type of analysis (this implementation recognizes only 'up'

if nargin<3
    mod='up';
end
nmc=[];
load (fl)

h=figure;
rt=Rate(find (EVENTS(1,:))/1000,10)/10;
plot (rt);
hold on;
[x1 ~]=ginput(1);
line([x1 x1],[min(rt) max(rt)]);
[x2 ~]=ginput(1);
close (h);
MX=round(10000*max(x1,x2));
MX=min (MX,length(EVENTS));
MN=round(10000*min(x1,x2));
MN=max(MN,1);
for i=1:4
    EVS{i}=find (EVENTS(i,MN:MX))+MN-1;
end

%%%%%%%%% CSP
ev=EVS{3}/1000;
rw=EVS{2}/1000;
lk=EVS{1}/1000;

[sel, ~]=listdlg('PromptString','Choose a method:','SelectionMode','single','ListString',{'PCA'; 'Value at +0.7V';...
    'Average +0.6V - +0.8V'; 'Normal distribution average'});


tr=repmat(TRI(1:2000,1),1,size(SWEEP_DATA,2));
CTR=SWEEP_DATA(1:2000,:)./tr;
CTR=mean(CTR(900:1500,:),1);

switch sel
    case 1
        %do PCA

    case 2
        dat=SWEEP_DATA(655,:);
    case 3
        dat=mean(SWEEP_DATA(640:680,:),1);
    case 4
        sh=pdf('norm',1:2000,660,15);
        s=size(SWEEP_DATA);
        for i=1:s(2)
            dat(i)=dot(SWEEP_DATA(1:2000,i),sh);
        end
end

rez=zeros(2,length(ev));
rez(:)=NaN;
rezc=zeros(1,pnts(2)*2);
nn=0;
for i=1:length(ev)
    f=find (SWEEP_TIMES>=ev(i));
    t1=SWEEP_TIMES(f(1));
    t2=SWEEP_TIMES(f(1)+pnts(2));
    if ~sum(lk>=t1 & lk<=t2)
        d=dat(f(1)+pnts(1):f(1)+pnts(2));
        dc=CTR(f(1)+pnts(1):f(1)+pnts(2));
        rez(2,i)=sum(lk>=t2 & lk<=t2+2)>0;
%         rezc(2,i)=rez(2,i);
        switch mod
            case 'up'
                rez(1,i)=(mean(d(1:2))-mean(d(end-1:end)))<0;
                ff=find(lk>=t2 & lk<t2+2);
                if ~isempty(ff)
                    ps=ff(1);
                    ff=find(SWEEP_TIMES>=lk(ps));
                    nn=nn+1;
                    rezc(nn,:)=dat(ff(1)-pnts(2)*2:ff(1)-1);
                end
%                 ft = fit((1:length(d))',d','poly1');
%                 a=coeffvalues (ft);
%                 rez(1,i)=a(1)>0;
%                 ft = fit((1:length(dc))',dc','poly1');
%                 a=coeffvalues (ft);
%                 rez(1,i)=a(1)>0;
        end
    end
end
nm(1)=sum(rez(2,find(rez(1,:)==1)))/length(find(rez(1,:)==1));
nm(2)=1-(length(find(rez(1,:)==0))-sum(rez(2,find(rez(1,:)==0))))/length(find(rez(1,:)==0));
% nmc(1)=sum(rezc(2,find(rezc(1,:)==1)))/length(find(rezc(1,:)==1));
% nmc(2)=(length(find(rezc(1,:)==0))-sum(rezc(2,find(rezc(1,:)==0))))/length(find(rezc(1,:)==0));
end