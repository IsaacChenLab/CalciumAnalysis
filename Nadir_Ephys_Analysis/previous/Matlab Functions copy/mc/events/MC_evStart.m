function MC_evStart(file,Time)

global CAT;
global DATE;
global DATA_DIR;

if nargin<2
    Time=nan;
end

[fileName,file]=MC_digFileName(file);
if exist(fileName,'file')
    fileName=[fileName ' -append'];
end

clear Ev*;
[BitTimes,BitMat,FileTimes]=MC_evGet(file,Time);

if issparse(BitMat)
    s=sprintf('save %s Bit* FileTimes;',fileName);
    eval(s);    
end

return;


