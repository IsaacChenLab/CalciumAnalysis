function DA_Prepare (fl)

% [a b]=butter(5,band*2/500);
% lfp=filtfilt(a,b,lfp);

fl=which(fl);
load (fl);
if ~exist('SWEEP_DATA','var') || ~exist('TRI','var')
    disp ('no variables saved');
    return
end
f=find(fl=='\');
ffl=[fl(1:3) 'BHV' fl(f(end):end)];
if fl(2)==':'
    nfl=[fl(1:end-4) '_prepared.mat'];
else
    nfl=['h:\Transfer\DA\DATA\' fl(1:end-4) '_prepared.mat'];
end
s=size(SWEEP_DATA);
c1=mean(SWEEP_DATA(4501,2:end))-mean(SWEEP_DATA(4000,2:end));
c1=c1/(TRI(4005,1)-TRI(5000,1));
% for i=2:s(2)
%     SWEEP_DATA(:,i)=SWEEP_DATA(:,i)./c1+TRI(1:5000,2);
% end
for i=2:s(2)
    c1=(SWEEP_DATA(4020,i)-SWEEP_DATA(4000,i))/(TRI(4005,1)-TRI(5000,1));
    SWEEP_DATA(:,i)=SWEEP_DATA(:,i)./c1+TRI(1:5000,2);
end
SWEEP_DATA=SWEEP_DATA(:,2:end);

save (nfl, 'SWEEP_DATA', 'SWEEP_TIMES', 'TRI');


a=questdlg('Do you need to synchronize behavioral events?','Synchronize','Yes','No','No');
if strcmp(a,'No')
    save (nfl,'EVENTS','-append');
    return;
end
load(ffl);
scrsz = get(0,'ScreenSize');
f=figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]);
c1=find(EVENTS)/1000;
ff=diff(c1);
c1=c1(ff>0.1);
c2=events{1};
% nnn=min(length(c1),length(c2));
nnn=10;
pos=1;
rep=1;
while rep
    hold off;
    plot (c1(1:nnn)-c1(pos),1,'kx');
    hold on;
    plot (c2(1:nnn)-c2(1),1,'ro');
    p=ginput(1);
    if isempty(p)
        rep=0;
    elseif p(1)<c2(5)-c2(1)
        pos=max(1,pos-1);
    else
        pos=min(length(c2),pos+1);
    end
end
close (f);
f=find(EVENTS);
ref1=f(pos);

if strcmp(Proc,'TL')
    evnts{1}=licks;
    evnts{2}=events{1};
    evnts{3}=events{2};
    evnts{4}=events{3};
    evnts{5}=events{10};
    evnts{6}=events{11};
    evnts{7}=events{12};
    evnts{8}=events{13};
    ref2=events{1}(1);
    EV=sparse(zeros(8,1));
    for i=1:8
       pos=ref1+round(1000*(evnts{i}-ref2));
       pos=pos(pos>0);
       EV(i,pos)=1;
    end
    EVT{1}=find(EV(1,:))./1000;
    EVT{2}=find(EV(2,:))./1000;
    EVT{5}=find(EV(5,:))./1000;
    EVT{6}=find(EV(6,:))./1000;
    EVT{7}=find(EV(7,:))./1000;
    EVT{8}=find(EV(8,:))./1000;
    ctp_l=XCorr (EVT{5},EVT{1},[-1 5],.1);
    ctm_l=XCorr (EVT{6},EVT{1},[-1 5],.1);
    ctp_r=XCorr (EVT{5},EVT{2},[-1 5],.1);
    clp_l=XCorr (EVT{7},EVT{1},[-1 5],.1);
    clm_l=XCorr (EVT{8},EVT{1},[-1 5],.1);
    clp_r=XCorr (EVT{7},EVT{2},[-1 5],.1);
    f=figure;
    sft=0;
    sfl=0;
    rep=1;
    while rep
        subplot (2,2,1); hold off;plot (0,0);
        rectangle ('Position', [10+sft 0 10 max(ctp_l)], 'FaceColor','y','EdgeColor','y');hold on;
        plot(ctp_l,'r'); plot (ctm_l,'k');
        subplot (2,2,2); hold off;plot (0,0);
        rectangle ('Position', [10+sfl 0 10 max(clp_l)], 'FaceColor','y','EdgeColor','y');hold on;
        plot(clp_l,'r'); plot (clm_l,'k');
        subplot (2,2,3); hold off;plot (0,0);
        rectangle ('Position', [10+sft 0 10 max(ctp_l)], 'FaceColor','y','EdgeColor','y');hold on;
        plot (ctp_r);
        subplot (2,2,4); hold off;plot (0,0);
        rectangle ('Position', [10+sfl 0 10 max(clp_l)], 'FaceColor','y','EdgeColor','y');hold on;
        plot (clp_r);
        [x y b]=ginput (1);
        if isempty(x)
            rep=0;
        else
            if b==1
                if x>10+sft && sft<100
                    sft=sft+1;
                elseif x<10+sft && sft>-9
                    sft=sft-1;
                end
            elseif b==3
                if x>10+sfl && sfl<100
                    sfl=sfl+1;
                elseif x<10+sfl && sfl>-9
                    sfl=sfl-1;
                end
            end
        end
    end
    close (f)
    EVT{5}=EVT{5}+sft/10;
    EVT{6}=EVT{6}+sft/10;
    EVT{7}=EVT{7}+sfl/10;
    EVT{8}=EVT{8}+sfl/10;
    EV(5,:)=0; EV(5,round(EVT{5}*1000))=1;
    EV(6,:)=0; EV(6,round(EVT{6}*1000))=1;
    EV(7,:)=0; EV(7,round(EVT{7}*1000))=1;
    EV(8,:)=0; EV(8,round(EVT{8}*1000))=1;
else
    evnts{1}=licks;
    evnts{2}=events{1};
    evnts{3}=events{2};
    evnts{4}=events{3};
    if ~isempty(events{9})
        evnts{5}=events{9};
    end
    ref2=events{1}(1);

    EV=sparse(zeros(4,1));
    for i=1:length(evnts)
       pos=ref1+round(1000*(evnts{i}-ref2));
       pos=pos(pos>0);
       EV(i,pos)=1;
    end
    
    EVT{1}=find(EV(1,:))./1000;
    EVT{2}=find(EV(2,:))./1000;
    EVT{3}=find(EV(3,:))./1000;
    EVT{4}=find(EV(4,:))./1000;
%     EVT{5}=find(EV(5,:))./1000;
    ctp_l=XCorr (EVT{3},EVT{1},[-1 5],.1);
    ctm_l=XCorr (EVT{4},EVT{1},[-1 5],.1);
    ctp_r=XCorr (EVT{3},EVT{2},[-1 5],.1);
    f=figure;
    sft=0;
    sfl=0;
    rep=1;
    while rep
        subplot (2,2,1); hold off;plot (0,0);
        rectangle ('Position', [10+sft 0 10 max(ctp_l)], 'FaceColor','y','EdgeColor','y');hold on;
        plot(ctp_l,'r'); plot (ctm_l,'k');
        subplot (2,2,3); hold off;plot (0,0);
        rectangle ('Position', [10+sft 0 10 max(ctp_l)], 'FaceColor','y','EdgeColor','y');hold on;
        plot (ctp_r);
        [x y b]=ginput (1);
        if isempty(x)
            rep=0;
        else
            if b==1
                if x>10+sft && sft<100
                    sft=sft+1;
                elseif x<10+sft && sft>-9
                    sft=sft-1;
                end
            elseif b==3
                if x>10+sfl && sfl<100
                    sfl=sfl+1;
                elseif x<10+sfl && sfl>-9
                    sfl=sfl-1;
                end
            end
        end
    end
    close (f)
    EVT{3}=EVT{3}+sft/10;
    EVT{4}=EVT{4}+sft/10;
    EV(3,:)=0; EV(3,round(EVT{3}*1000))=1;
    EV(4,:)=0; EV(4,round(EVT{4}*1000))=1;
%     EV(5,:)=0; EV(5,round(EVT{5}*1000))=1;
    
end
EVENTS=sparse(EV);
save (nfl,'EVENTS','-append');