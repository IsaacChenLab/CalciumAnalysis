function DIY_Ca_Event_Inspector(ROI,num_frames,fps)

fprintf("\tSelect folder containing data to be analyzed...");
data_path = uigetdir('','Select folder containing data to be analyzed');
fprintf("Selected!\n");

fprintf("\tSelect folder to save figures to...");
save_path = uigetdir('','Select folder to save figures to');
fprintf("Selected!\n");

%load data file
load(strcat(data_path, '/processed_analysis.mat'));     

%collect time axis
time = 1:num_frames;
time = time ./ fps;

for k = 1:length(ROI)

    %collect deltaf trace and raw trace
    dF_trace = processed_analysis.dF_cell(ROI(k),:);
    raw_trace = processed_analysis.F_cell(ROI(k),:);

    %collect spike times
    spike_times = processed_analysis.Spikes_cell{ROI(k)};
    spike_times = spike_times ./ fps;

    %create the figure and axes
    name = ['ROI ' num2str(ROI(k))];
    f = figure('Name', name, 'NumberTitle','off');
    ax1 = axes('Position',[0.1 0.6 0.8 0.3]);
    ax2 = axes('Position',[0.1 0.1 0.8 0.3]);

    % plot each trace
    plot(ax1,time,raw_trace,'k');
    plot(ax2,time,dF_trace,'k');

    %label axes
    xlabel(ax1,'Time (s)');
    ylabel(ax1,'Calcium fluorescence (a.u.)');
    title(ax1, name);

    xlabel(ax2,'Time (s)');
    ylabel(ax2,'\DeltaF/F');

    %add blue lines for each event to both axes
    axes(ax1);
    for i=1:length(spike_times)
        line([spike_times(i),spike_times(i)], ylim(ax1));
    end
    
    axes(ax2);
    for i=1:length(spike_times)
        line([spike_times(i),spike_times(i)], ylim(ax2));
    end

    %save figure as a jpeg
    saveas(f,strcat(save_path, '/', name, ' - 85', '.jpg'));
    
end
fprintf("\tAll done - files are saved\n");

end