function BHV_Compact (nm)

al=dir([nm '*']);
if isempty(al)
    error ('no such name');
end

dt=datenum(clock);
el=zeros(length(al),1);
for i=1:length(al)
    el(i)=dt-al(i).datenum;
end
[~, pos]=sort(el);
pos=pos(end:-1:1);
v=[];
for i=1:length(al)
    Proc=[];
    load (al(pos(i)).name,'Proc');
    v{i}=([ num2str(i) '     Date: ' al(pos(i)).date '; Procedure: ' Proc]);
end
q=Confirm_dlg(v);
if strcmp(q,'No')
    error ('please fix the file names / dates');
end

vv=inputdlg({'Combine sessions (eg. 1:5) #:'; 'Licking period (sec):'},'Parameters',2,{'1:5'; '200'});
val=str2double(vv{2});
f=find(vv{1}==':');
strt=str2double(vv{1}(1:f-1));
if strcmp(vv{1}(f+1:end),'end')
    strt(2)=length(al);
else
    strt(2)=str2double(vv{1}(f+1:end));
end
if isnan(val)
    val=200;
end
EVNTS=[];
EVNTS{100}=[];
LCKS=[];
ref=0;
for ii=strt(1):strt(2)
    i=[];
    load (al(pos(ii)).name);
    events=events(i,:);
    LCKS=[LCKS licks+ref];
    licks(end+1)=5000;
    df=diff(licks);
    fdf=find(df>val);
    outs=licks(fdf+1);
    outs(2,:)=licks(fdf);
    
    for z=1:length(events)
        for y=1:length(events{z})
            s=(outs(2,:)<events{z}(y)) & (outs(1,:)>events{z}(y));
            if ~sum(s(:))
                EVNTS{z}(end+1)=events{z}(y)+ref;
            end
        end
    end
    ref=ref+5000;
end
nnm=[nm '_Sessions_' num2str(strt(1)) '-' num2str(strt(2)) '.mat'];
save (nnm,'EVNTS', 'LCKS');