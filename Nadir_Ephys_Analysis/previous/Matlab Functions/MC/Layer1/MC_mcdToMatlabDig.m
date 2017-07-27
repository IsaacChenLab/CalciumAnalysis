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
    flname=sprintf('%s\\PreparedFiles\\%s_%s_%d',DATA_DIR,CAT,DATE,file(i));
    %if exist(flname,'file')
        flnamec=[flname ' -append'];
    %end
    % digital
    [digData,xtime]=MC_getDigital(fileName);
    s=sprintf('save %s digData;',flnamec);
    eval(s);
end
close(h1);

return;
