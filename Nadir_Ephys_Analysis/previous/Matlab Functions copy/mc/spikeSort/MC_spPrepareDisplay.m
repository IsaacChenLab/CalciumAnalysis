function MC_spPrepareDisplay(file,chs)


global CAT;
global DATE;
global DATA_DIR;
FROM=20e3; TO=60e3;


d=1;
for ch=chs
    d
    clear dispData*;
    dispData=[]; xtime=[];
    for j=1:length(file)
        j
        fileName=MC_fileName(file(j));
        [dataStream,totalMS,Hz]=MC_openFile(fileName);
        if totalMS>TO,
            [rd,xt]=MC_getElectrodes(fileName,ch,[FROM TO]);
            dispData=[dispData rd];
        end
    end
    s=sprintf('dispData_%d=dispData;',ch);
    eval(s);
    cd(DATA_DIR);
    cd('DisplayFiles');
    if d==1,
        s=sprintf('save D%s_%s_%d-%d dispData_%d;',CAT,DATE,file(1),file(end),ch);
        eval(s);
    else
        s=sprintf('save D%s_%s_%d-%d -append dispData_%d;',CAT,DATE,file(1),file(end),ch);
        eval(s);
    end

    d=d+1;
end
cd(DATA_DIR);
cd('DisplayFiles');
s=sprintf('save D%s_%s_%d-%d -append xt Hz;',CAT,DATE,file(1),file(end));
eval(s);


return;
