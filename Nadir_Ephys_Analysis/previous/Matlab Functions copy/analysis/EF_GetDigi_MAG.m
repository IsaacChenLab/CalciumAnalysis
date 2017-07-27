function EF_GetDigi_MAG (fl)

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
DIG(round(1000*xs(df<0)))=1;
ff=find(DIG);
df=diff(ff);
DIG_UV=sparse(zeros(1,round(1000*length(data)/freq)));
DIG_UV(ff(find (df<5000)+2))=1;
DIG_LIGHT=DIG;
DIG_LIGHT(ff(find (df<5000)+1))=0;
DIG_GREEN=sparse(zeros(1,round(1000*length(data)/freq)));
DIG_GREEN(ff(find (df<5000)-1))=1;


save (fl, 'DIG', 'DIG_UV', 'DIG_LIGHT', 'DIG_GREEN', '-append');