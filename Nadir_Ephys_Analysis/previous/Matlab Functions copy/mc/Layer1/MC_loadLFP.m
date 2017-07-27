function [lfp,totalMs,LFP_Hz,LFP_time]=MC_loadLFP(file,ch)

global CAT;
global DATE;
global DATA_DIR;

cd(DATA_DIR);
cd('lfpFiles2');

s=sprintf('load L2%s_%s_%d lfpData_%d totalMs LFP_Hz;',CAT,DATE,file,ch);
eval(s);


s=sprintf('lfp=lfpData_%d; clear lfpData_%d;',ch,ch);
eval(s);

if nargout==4
    s=sprintf('load L2%s_%s_%d LFP_time;',CAT,DATE,file);
    eval(s);
end

return;





