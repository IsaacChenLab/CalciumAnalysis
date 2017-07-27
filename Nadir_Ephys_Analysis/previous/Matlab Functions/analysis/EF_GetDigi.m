function EF_GetDigi (fl)

try
    load(fl)
    freq=info.Fs;
    data=digi;
    clear info digi;
catch
    warndlg('Incorrect or missing file');
    return;
end

xs=0:1/freq:(length(data)-1)/freq;
df=diff(data);
DIG=sparse(zeros(1,round(1000*length(data)/freq)));
DIGe=sparse(zeros(1,round(1000*length(data)/freq)));
DIG(round(1000*xs(df<0)))=1;
DIGe(round(1000*xs(df>0)))=1;
ff=find(DIG);
df=diff(ff);
DIGs=sparse(zeros(1,round(1000*length(data)/freq)));
DIGs(ff([1 find(df>1000)+1]))=1;
fs=find(DIGs);
fe=find(DIGe);
ffe=find(diff(fe)>1000);
PT=DIG(fs(1):fe(ffe(1)));
PATT=0.001:1/1000:length(PT)/1000;
PATT(2,:)=PT;
PATT=full(PATT);
q=questdlg ('Synchronize behavior?','Question','Yes','No','No');
if strcmp(q,'Yes')
    fl=which(fl);
    f=find(fl=='\');
    nfl=[fl(1:f(1)) 'BHV\' fl(f(3)+1:f(4)-1) '.mat'];
    try
        load (nfl);
        EVENTS=[];
        EVENTS{20}=[];
        pdig=find(DIGs);
        plas=events{9}*1000; %%%
        rep=1;
        f=figure;
        df=0;
        while rep
            hold off
            plot (pdig(1:10+df)-pdig(1+df),0,'ro');
            hold on;
            plot (plas(1:10)-plas(1),0,'xk');
            [x y bt]=ginput(1);
            if isempty(x)
                rep=0;
            else
                if bt==1 && df>0
                    df=df-1;
                elseif bt==3 && df<20
                    df=df+1;
                end
            end
        end
        close (f);
        mn=events{9}(1); %%%
        pl=pdig(1+df);
        COMMENTS={'LICKS', 'REWARD', 'CSP', 'CSM', 'LASER0', 'LASER1','LASER_ALL'};
        if isempty(licks)
            EVENTS{1}=[];
        else
            EVENTS{1}=(licks-mn)*1000+pl;
        end
        if isempty(events{1})
            EVENTS{2}=[];
        else
            EVENTS{2}=(events{1}-mn)*1000+pl;
        end
        if isempty(events{2})
            EVENTS{3}=[];
        else
            EVENTS{3}=(events{2}-mn)*1000+pl;
        end
        if isempty(events{3})
            EVENTS{4}=[];
        else
            EVENTS{4}=(events{3}-mn)*1000+pl;
        end
        if isempty(events{10})
            EVENTS{5}=pdig;
        else
            EVENTS{5}=pdig(~events{10});
            EVENTS{6}=pdig(events{10}==1);
        end
        EVENTS{7}=find (DIG);
    catch
        error (['No behavior file named ' nfl])
    end
    save (fl, 'DIG', 'DIGe', 'DIGs', 'EVENTS', 'COMMENTS', 'PATT', '-append');
end
save (fl, 'DIG', 'DIGe', 'DIGs', 'PATT', '-append');