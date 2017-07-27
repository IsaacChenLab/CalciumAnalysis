function EF_AverageRaws (list)

if nargin<1
    list=1:16;
end
fl=dir('raws*');
if isempty(fl)
    return;
end
load (['raws' num2str(list(1),'%02.f')]);
raw0=raw;
h=waitbar(0,'...');
for i=2:length(list)
    waitbar(i/length(fl));
    load (['raws' num2str(list(i),'%02.f')]);
    raw0=raw0+raw;
end
close (h);
raw=raw0./length(fl);
save('raws00.mat','raw','info');
