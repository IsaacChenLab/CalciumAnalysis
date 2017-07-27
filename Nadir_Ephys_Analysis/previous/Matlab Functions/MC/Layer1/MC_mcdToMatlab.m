function MC_mcdToMatlab(file,chs)

global CAT;
global DATE;
global DATA_DIR;

fileName=MC_fileName(file);
file=MC_fileNumber(fileName); 
h1=waitbar(0,'Preparing files');
for i=1:length(file)
    fileName=MC_fileName(file(i));
    [dataStream,totalMs,Hz]=MC_openFile(fileName);
    flname=sprintf('%sPreparedFiles\%s_%s_%d',DATA_DIR,CAT,DATE,file(i));
    if exist(flname,'file')
        flnamec=[flname ' -append'];
    end
    % digital
    [digData,xtime]=MC_getDigital(fileName);
    s=sprintf('save %s digData;',flnamec);
    eval(s);
    flnamec=[flname ' -append'];
    d=1;
    for ch=chs
        waitbar(((i-1)*length(chs)+d)/(length(chs)*length(file)),h1);
    % get raw data
        clear rawData* xtime;
        s=sprintf('[rawData_%d,xtime]=MC_getElectrodes(fileName,ch);',ch);
        eval(s);
        s=sprintf('save %s rawData_%d;',flnamec,ch);
        eval(s);
        d=d+1;
    end
    % save time
    s=sprintf('save %s xtime totalMs Hz;',flname);
    eval(s);
end
close(h1);

return;
