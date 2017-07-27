function MC_scriptSharpWave(file,ch,Time)

global CAT;
cdc;
if nargin<3
    Time=nan;
end

sharp_file=sprintf('Sharp_%s.mat',CAT);
if exist(sharp_file,'file')
    eval(['load ' sharp_file]);
end

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);    
nfl=sprintf('%d',file);

s=sprintf('[SharpTimes_%s_%d,SharpMat_%s_%d]=MC_findSharpWaves(file,ch,Time);',nfl,ch,nfl,ch);
eval(s);

s=sprintf('save %s Sh*;',sharp_file);
eval(s);    

return;
