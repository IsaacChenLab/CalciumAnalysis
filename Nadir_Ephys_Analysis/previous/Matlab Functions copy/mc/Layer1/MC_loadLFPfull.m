function [lfp,totalMs,LFP_Hz]=MC_loadLFPfull(file,ch)

global CAT;
global DATE;
global DATA_DIR;

cd(DATA_DIR);
cd('lfpFiles2');

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  
lfp=[];
for fl=file
    s=sprintf('load L2%s_%s_%d lfpData_%d totalMs LFP_Hz;',CAT,DATE,fl,ch);
    eval(s);
    s=sprintf('lfp=[lfp; lfpData_%d]; clear lfpData_%d;',ch,ch);
    eval(s);
end

return;






