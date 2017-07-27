function MC_spPrepareFilesDig(file)

global CAT;
global DATE;
global DATA_DIR;

fileName=MC_fileName(file);
file=MC_fileNumber(fileName); 
h1=waitbar(0,'Preparing files');
for i=1:length(file)
    fileName=MC_fileName(file(i));
    [dataStream,totalMs,Hz]=MC_openFile(fileName);
    cd(DATA_DIR);
    cd('PreparedFiles');
    % dig data
    [digData]=MC_getDigital(fileName);
    s=sprintf('save %s_%s_%d -append digData;',CAT,DATE,file(i));
    eval(s);
end
close(h1);

return;
