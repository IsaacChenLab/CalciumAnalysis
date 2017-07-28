function Waveform_Library_Editer()

% Used for adding waveforms to library used in event detection

fprintf("\nSelect spikes.mat file...");
[data_file, data_path] = uigetfile('*.mat','Select spikes.mat file...');
fprintf("Selected!\n");
spikesFile = [data_path, data_file];
load(spikesFile);

fprintf("Select processed_analysis.mat file...");
[data_file, data_path] = uigetfile('*.mat','Select processed_analysis.mat file...');
fprintf("Selected!\n");
load([data_path, data_file]);

ROIs = input('Provide vector of ROIs to inspect: ');
fprintf("\n");

for r = ROIs
    
    h = figure;
    plot(processed_analysis.dF_cell(r,:),'k');
    xlabel('Frame #');
    ylabel('DeltaF/F');
    title(sprintf('ROI %.0f\nSelect two points that mark the start and end of a waveform to add to the library\nPress <Enter> after selecting the second point', r));
    
    useThisOne = input(['ROI ' num2str(r) ': use = 1, next ROI = 0: ']);
    while useThisOne
        plot(processed_analysis.dF_cell(r,:),'b');
        [x,~] = getpts(h);
        if length(x) == 2
            spikes{end+1} = processed_analysis.F_cell(r,floor(x(1)):floor(x(2)));
            plot(processed_analysis.dF_cell(r,:),'k');
            fprintf("\t---Wave #%.0f added---\n\t",length(spikes));
        else
            fprintf("\tBad attempt\n\t");
        end
        useThisOne = input('another = 1, next ROI = 0: ');
    end
    
    close(h);
end

% display library with all waves added
fprintf("\nDisplaying current library for inspection\n");
f = PlotLibrary(spikes);

% remove undesired waves
badWaves = input('List in a vector which waves should be removed ([] = none) ');
for b = length(spikes):-1:1
    if sum(badWaves == b) > 0
        spikes(:,b) = [];
    end
end

% display final library and save
close(f);
fprintf("Displaying final library\n");
PlotLibrary(spikes);
save(spikesFile,'spikes');

end


% helper function for plotting the library
function f = PlotLibrary(spikes)
n = length(spikes);
r = ceil(sqrt(n));
c = ceil(n/r);
f = figure;
for i = 1:n
    subplot(r,c,i);
    plot(spikes{i},'k');
    title(['Waveform #' num2str(i)]);
end
end
