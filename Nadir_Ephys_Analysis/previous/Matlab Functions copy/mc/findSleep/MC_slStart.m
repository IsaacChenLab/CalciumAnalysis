function MC_slStart(file,ch,eog_ch,addrem)

global CAT;
global DATE;
global DATA_DIR;

if nargin<3
    eog_ch=nan;
end
if nargin<4
    addrem=0;
end

[fileName,file]=MC_digFileName(file);
if exist(fileName,'file')
    fileName=[fileName ' -append'];
end

[swsMat,swsTimes,remMat,remTimes,FileTimes]=MC_slFindSleep(file,ch,eog_ch);

if (~isnan(swsMat) | ~isnan(remMat)) & (~isempty(swsTimes) | ~isempty(remTimes))
    if ~addrem,
        s=sprintf('save %s swsMat swsTimes remMat remTimes FileTimes;',fileName);
    else
        s=sprintf('save %s remMat remTimes;',fileName);        
    end
    eval(s);    
end

return;


