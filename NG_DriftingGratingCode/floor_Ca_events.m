%changes the data file so that the timing for each event is a whole number

load('processed_analysis.mat');

for i = 1:length(processed_analysis.Spikes_cell)
    processed_analysis.Spikes_cell{i} = floor(processed_analysis.Spikes_cell{i});    
end

save('processed_analysis.mat', '-append');