load('/Users/Nick/Documents/ChenLab/ChenLabData/BCT debugging/test/processed_analysis.mat')

spikes = processed_analysis.Spikes_cell;
event_count = sum(processed_analysis.dat(:,2));

x = zeros(event_count,1);
y = zeros(event_count,1);

p = 1;
for i = 1:length(spikes)
    for j = 1:length(spikes{i})
        x(p) = spikes{i}(j);
        y(p) = i;
        p = p + 1;
    end
end

x = x ./ 5;

f = figure; %('Visible', 'off');
scatter(x,y,20,'filled');
  
xlabel('Time (s)');
ylabel('Cell number');
title('Ca Event Scatter plot');

%saveas(f,'/Users/Nick/Desktop/histo1.jpg');