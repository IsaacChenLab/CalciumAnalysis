function MC_slStart(file,ch)

if nargin<2
    error('not enough args');
end
global CAT;
global DATE;
global DATA_DIR;

[fileName,file]=MC_digFileName(file);
if exist(fileName,'file')
    fileName=[fileName ' -append'];
end

    clear sh*;
    [sharpMat,sharpTimes,FileTimes]=MC_swFindSharp(file,ch);
    s=sprintf('sharpMat_%d=sharpMat; sharpTimes_%d=sharpTimes;',ch,ch);
    eval(s);
    if issparse(sharpMat)
        s=sprintf('save %s sharpMat_%d sharpTimes_%d FileTimes;',fileName,ch,ch);
        eval(s);    
    end

return;


