function MC_spPrepareFiles(file,chs)

global CAT;
global DATE;
global DATA_DIR;

fileName=MC_fileName(file);
file=MC_fileNumber(fileName); 
h1=waitbar(0,'Preparing files');
for i=1:length(file)
    fileName=MC_fileName(file(i));
    [dataStream,totalMs,Hz]=MC_openFile(fileName);
    d=1;
    flfs=sprintf('%s_%s_%d.mat',CAT,DATE,file(i));
    for ch=chs
        waitbar(((i-1)*length(chs)+d)/(length(chs)*length(file)),h1);
    % get raw data
        clear rawData* xtime;
        s=sprintf('[rawData_%d,xtime]=MC_getElectrodes(fileName,ch);',ch);
        eval(s);
        cd(DATA_DIR);
        cd('PreparedFiles');
        if ~exist(flfs,'file'),
            s=sprintf('save %s rawData_%d;',flfs,ch);
            eval(s);
        else
            s=sprintf('save %s -append rawData_%d;',flfs,ch);
            eval(s);
        end
        d=d+1;
    end
    % dig data
    [digData]=MC_getDigital(fileName);
    s=sprintf('save %s -append digData;',flfs);
    eval(s);
    % save time
    s=sprintf('save %s -append xtime totalMs Hz;',flfs);
    eval(s);
end
close(h1);

return;
