load('/Users/Nick/Documents/ChenLab/ChenLabData/BCT debugging/test/processed_analysis.mat')

events = processed_analysis.dat(:,2);

%edges = [0,1,4,7,10,max(events)];
f = figure; %('Visible', 'off');
histogram(events,max(events));

xlabel('Events per cell');
ylabel('Number of cells');
title('Event Histogram');

%saveas(f,'/Users/Nick/Desktop/histo1.jpg');