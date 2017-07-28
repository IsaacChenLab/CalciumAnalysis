load('/Users/Nick/Documents/ChenLab/ChenLabData/BCT debugging/test/processed_analysis.mat')

spikes = processed_analysis.Spikes_cell;

dlmwrite('/Users/Nick/Desktop/all_events.csv',[]);

fid = fopen(['/Users/Nick/Desktop/all_events.csv'],'w');
fprintf(fid,'ROI, Event timings (s)\n');
fclose(fid);
        
for x = 1:length(spikes)
    labelled_spikes = [x spikes{x}./5];
    dlmwrite('/Users/Nick/Desktop/all_events.csv', labelled_spikes, '-append');
end
