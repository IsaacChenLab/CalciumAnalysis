function MC_lfpPrepareFiles(file,chs)

global CAT;
global DATE;
global DATA_DIR;

R1=25;
R2=50;
disp(DATE);
for i=1:length(file)
    d=1;
    [rawData,xtime,Hz,totalMs]=MC_loadElectrodes(file(i),1);
        
    LFP_time=decimate(xtime,R1);
    LFP_Hz=Hz/R1;
    % save time
    s=sprintf('save L%s_%s_%d -append LFP_time totalMs LFP_Hz;',CAT,DATE,file(i));
    cd(DATA_DIR);
    cd('lfpFiles');
    eval(s);
    
    LFP_time=decimate(xtime,R2);
    LFP_Hz=Hz/R2;
    s=sprintf('save L2%s_%s_%d -append LFP_time totalMs LFP_Hz;',CAT,DATE,file(i));
    cd(DATA_DIR);
    cd('lfpFiles2');
    eval(s);
end

return;
