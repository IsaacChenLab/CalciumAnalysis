function rez=temp_PlotAllSpikes (dr, evs, wnd, bin)

cd(['Z:\EPhys\To Matlab\' dr]);
fd=dir('*_spk.mat');
figure;
col={'k' 'r' 'b' 'g'};
rez=[];
for i=1:length(fd)
    load (fd(i).name);
    tf=find(fd(i).name=='_');
    nr=str2double(fd(i).name(tf(end)-2:tf(end)-1));
    subplot (4,4,nr);
    hold on
    for j=1:size(SpikeTimes,2)
        plot (EF_XCorr(evs, SpikeTimes{j}, wnd, bin, 1000),col{j});
        rez{end+1}=EF_XCorr(evs, SpikeTimes{j}, wnd, bin, 1000);
    end
    title (['Ch. ' num2str(nr)]);
end