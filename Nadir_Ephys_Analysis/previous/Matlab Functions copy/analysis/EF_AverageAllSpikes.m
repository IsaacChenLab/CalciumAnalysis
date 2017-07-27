function [mat] = EF_AverageAllSpikes (ref, wnd, bin,freq)

mat=[];
d=dir('*_spk.mat');
if isempty(d)
    return;
end
load (d(1).name,'SpikeTimes');
mat=EF_XCorr (ref,SpikeTimes{1},wnd,bin,freq);
for i=1:length(d)
    load (d(i).name,'SpikeTimes');
    for j=1:size(SpikeTimes,2)
        mat(:,end+1)=EF_XCorr (ref,SpikeTimes{j},wnd,bin,freq);
    end
end
mat=mat(:,2:end);