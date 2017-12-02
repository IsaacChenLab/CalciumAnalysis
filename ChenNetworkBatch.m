function selectedFolders = ChenNetworkBatch(outputFolder, real_FPS, FC_method, FC_inactive, FC_islands)
%Computes functional network analysis metrics

%INPUT
% outputFolder - name of the folder in which all output will be deposited
% real_FPS - video frame rate, must provide this parameter
% FC_method - set to 'raw' to run raw trace method for FC -- otherwise
%     FLUROSNNAP default method will be used
% FC_inactive - set to 'include inactive' to include cells without events in
%     FC analysis, otherwise those cells will be ignored
% FC_islands - set to 'remove islands' to remove from network any cell that
%      has no edges

%set(0,'DefaultFigureVisible','off');

if ~exist('real_FPS', 'var')
    real_FPS = 5;
end
load('params.mat');
params.fps = real_FPS;


% get the spike library
fprintf('\nSelect spike library to be used...');
[spike_file, spike_path] = uigetfile('*.mat', 'Select .mat file');
fprintf([spike_file ' was selected!\n\n']);
params.spikeLib = fullfile(spike_path, spike_file);


%set parameters appropriately
if exist('FC_method', 'var') ~= 0
    params.FC_method = FC_method;
else
    params.FC_method = 'FluoroSNNAP default';
end

if exist('FC_inactive', 'var') ~= 0
    params.FC_inactive = FC_inactive;
else
    params.FC_inactive = 'exclude inactive';
end


if exist('FC_islands', 'var') == 0
    FC_islands = 'keep islands';
elseif ~strcmpi(FC_islands, 'remove islands')
    FC_islands = 'keep islands';
end

%feedback on the parameters
fprintf('\tFrames per second = %f\n', real_FPS);
if strcmpi(params.FC_method, 'raw')
    fprintf('\tFC_method = raw trace\n');
    if strcmpi(params.FC_inactive, 'include inactive')
        fprintf('\tFC_inactive = INCLUDING cells that have zero Ca events in FC analysis\n');
    else
        fprintf('\tFC_inactive = EXCLUDING cells that have zero Ca events in FC analysis\n');
    end
else
    fprintf('\tFC_method = FluoroSNNAP default\n');
    fprintf('\tFC_inactive = EXCLUDING cells w/o Ca events (required by the deafault FC method)\n');
end

if strcmpi(FC_islands, 'remove islands')
    fprintf('\tFC_islands = REMOVE all islands from the network\n');
else
    fprintf('\tFC_islands = KEEP islands in the network\n');
end

fprintf('\tNo threshold for weighted FC matrix being applied\n');
fprintf('\tCa event threshold = %.2f\n', params.event_thresh);

% turn FluoroSNNAP analysis modules on and off appropriately
params.analyze.deltaF = 1;
params.analyze.detect_events = 1;
params.analyze.FC = 1;
params.analyze.sca = 1;
params.analyze.controllability = 0;
params.analyze.kinetics = 0;
params.analyze.spike_probability = 0;
params.analyze.ensembles = 0;
params.analyze.figure = 0;
save('params.mat', 'params', '-append');


selectedFolders = {};


% GUI big folders full of videos all at once
getAnotherBig = input('\nDo you want to select a BIG folder containing lots of videos that all need to be analyzed? 0 = no, 1 = yes ');

while(getAnotherBig)
    bigFolder = uigetdir('', 'Select BIG folder containing small folders with tiff stack (.tif) and Seg file (.mat)');
    listing = dir(bigFolder)
    
    for ff = 3:length(listing)
        ff = listing(ff);
        videoFolder = fullfile(bigFolder, ff.name);
        if ff.isdir && ~isempty(dir(fullfile(videoFolder, '*.tif'))) && ~isempty(dir(fullfile(videoFolder, '*.mat')))
            selectedFolders{end+1} = videoFolder;
        end
    end
       
    getAnotherBig = input('Select another BIG folder? 0 = no, 1 = yes ');
end


% GUI for adding individual videos

getAnotherSmall = input('\nDo you want to add a small containing a single video? 0 = no, 1 = yes ');

while(getAnotherSmall)
    smallFolder = uigetdir('', 'Select small folder containing single video, ie tiff stack (.tif) and Seg file (.mat)');
  
    if ~isempty(dir(fullfile(smallFolder, '*.tif'))) && ~isempty(dir(fullfile(smallFolder, '*.mat')))
        selectedFolders{end+1} = smallFolder;
    end  

    if isempty(dir(fullfile(smallFolder, '*.tif')))
        disp (['There is no .tif file in the selected directory: ' smallFolder]);
    end
    if isempty(dir(fullfile(smallFolder, '*.mat')))
        disp (['There is no seg file (.mat) in the selected directory: ' smallFolder]);
    end
    
    getAnotherSmall = input('Select another single video? 0 = no, 1 = yes ');
end



% loop through files & analyze
% disp(['These folders have been selected:' strjoin(selectedFolders, '\n')]);
for zz = 1:length(selectedFolders)
    
    selected_folder = selectedFolders{zz};
    outputDir = fullfile(selected_folder, outputFolder);
    mkdir(outputDir);
    params.outputDir = outputDir;
    save('params.mat', 'params', '-append');
    
    disp(['Batch analysis begun on ' selected_folder]);
    
    % Run FluoroSNNAP ROI Analysis
    addpath('FluoroSNNAP_code');
    addpath('textprogressbar');
    addpath(fullfile('FluoroSNNAP_code', 'oopsi-master'));
    rawData = AnalyzeData_2(selected_folder, real_FPS);
    
    % Run FluoroSNNAP Processing if video has at least 10 frames
    if size(rawData{1}.F_cell,2) < 10
        continue;
    end
    
    processed_analysis = PostProcess_test(selected_folder);
    disp('Running BCT network analysis');
    
    % Extract FC matrix from FluoroSnnap

    load(fullfile(outputDir, 'processed_analysis.mat'));
    fc_matrix = processed_analysis.FC.CC.C;
    spikes = processed_analysis.Spikes_cell;
    
    % Check if video has any spikes, skip if not
    anySpikes = 0;
    for aa = 1:length(spikes)
        if ~isempty(spikes{aa})
            anySpikes = 1;
            break;
        end
    end
    if ~anySpikes
        continue;
    end
    
    % remove islands if called for
    if (strcmpi(FC_islands, 'remove islands'))
        fc_matrix( all(~fc_matrix,2), : ) = [];
        fc_matrix( :, all(~fc_matrix,1) ) = [];
    end
    
    % add BCT to MATLAB path
    addpath('BCT_code');
    
    % initialize output struct
    output = struct();
    
    % SYNCHRONIZATOIN
    % Not a network measure, but a measure of whether groups of neurons in
    % the organoid all fire together. GSI is between -1 and 1. Where 0 is 
    % total independence of all neurons, 1 is all neurons only fired
    % together. Negative numbers are anti-correlation and are unlikely to
    % be observed. Is computed by FluoroSNNAP instantaneous phase method,
    % not by BCT. Default is that a cluster must contain at least 3 ROI.
    output.global_synch_index = processed_analysis.SI;
    output.number_synch_clusters = processed_analysis.SynchroCluster.Num;
    
    % DENSITY
    % Density = fraction of present connections to possible connections.
    % # present connections = # nonzero elements in top triangle of matrix
    % # possible connections = ((# nodes) choose 2)
    [output.density,output.vertices,output.edges] = density_und(fc_matrix);
    
    % DEGREE DISTRIBUTION
    % Degree per cell is the number of connections it has
    % within the network Simply, this is the sum of pair-wise
    % connections (FC matrix elements) between Cell 1 and all other
    % cells that are not zero. This would be sensitive to a threshold
    % we could set in generating the FC matrix. Strength is the same as
    % degree but accounts for connection weights
    output.degree_distribution = (degrees_und(fc_matrix))';
    output.degree_avg = mean(output.degree_distribution);
    output.degree_avg_normalized = output.degree_avg/size(fc_matrix,1);
    output.strength_distribution = (strengths_und(fc_matrix))';
    output.strength_mean = mean(output.strength_distribution);
    
    % CLUSTERING COEFFICIENT
    % Clustering coeff of node n is the number of complete
    % triangles there are centered on nweighted by the geometric mean of
    % the edges involved. A complete triangle centered on n is when two
    % of n's neighbors are neighbors with each other. Network average is
    % the average of each node's clustering coeff. The transitivity is
    % the sum of the geometric mean of all triangles in the whole
    % network divided by the total number of potential triangles (pairs
    % of mutual friends) which is the sum of (degree(n) * degree(n-1))
    % summed over every node n
    output.clustering_coeff_distribution = clustering_coef_wu(fc_matrix);
    output.network_avg_clustering_coeff = mean(output.clustering_coeff_distribution);
    output.transitivity = transitivity_wu(fc_matrix);
    %[C_pos,C_neg,output.average_clustering_coeff,Ctot_neg] = clustering_coef_wu_sign(fc_matrix,1);
    
    % MODULARITY
    % Modularity of graph is the degree to which that graph
    % can be divided into a set of communities (modules) such that
    % within-community edges are maximized and between-community edges
    % are minimized. Finding the the modularity of the graph requires
    % finding the optimal community structure (ie the one that maximizes
    % modularity) which has been proven to be NP-hard. So this algorithm
    % provides a best estimate of the true modularity. Cj is the optimal
    % community structure (ie the actual modules).
    [Cj, output.modularity] = modularity_und(fc_matrix);
    
    
    % LENGTH MATRIX
    % Must convert the FC matrix into a 'length' matrix.
    % A length between two nodes is inverse of the connectivity weight.
    % A shorter length means two nodes are better connected.
    length_mat = weight_conversion(fc_matrix,'lengths');
    
    % DISTANCE MATRIX:
    % Each matrix entry (i,j) is the length of the shortest path
    % between nodes i and j NOTE: Previously the length matrix (where
    % edge weights are inversed so as to represent lengths) was
    % mistaken for the distance matrix. Bug fixed on 7/10/17 by NG
    dist_mat = distance_wei(length_mat);
    
    % GLOBAL EFFICIENCY
        % Average of inverse of distance between all pairs of
        % nodes in the network (fxnal segregation measure)
        % output.global_efficiency = mean(1./D);
    % CHARACTERISTIC PATH
        % Average distance between all pairs of nodes in
        % network (fxnal segregation measure)
    % NODE ECCENTRICITY
        % Eccentricity of node n is distance from n to the node that is
        % farthest from n
    % RADIUS
        % the minimum node eccentricity in the network
    % DIAMETER
        % the maximum node eccentricity in the network
    [output.characteristic_path,output.global_efficiency,output.eccentricities,...
        output.radius,output.diameter] = charpath(dist_mat);
    output.mean_eccentricity = mean(output.eccentricities);
    
    %ASSORTATIVITY
    % correlation coefficient for the strengths of nodes that are connected
    output.assortativity = assortativity_wei(fc_matrix, 0);
    %^^ 0 flag = undirected graph
    
    % RICH CLUB CURVE
    % For all possible degrees k, the fraction of edges that exist
    % between all nodes of degree k out of the maximum number of such
    % edges that could possibly exist returns a vector of length k
    % where k is the maximum degree of any node in the network
    output.rich_club_curve = (rich_club_wu(fc_matrix))';
    
    % BETWEENNESS CENTRALITY
    % Fraction of paths in network that contain each node.
    % The nodes with greatest centrality are 'hub' nodes. The algorithm
    % finds the number of all shortest paths between every set of two
    % nodes that passes through any one node. The greater number of
    % paths through the node, the greater centrality
    B = betweenness_wei(length_mat);
    %must normalize
    n = length(fc_matrix);
    output.betweenness_distribution = B./((n-1)*(n-2));
    
    %SAVE THE OUTPUT;
    %make the directory

    BCT_directory = fullfile(outputDir, 'BCT');
    if ~exist(BCT_directory, 'dir')
        mkdir(BCT_directory);
    end
    cprintf('*blue','%s\n', ['Saving BCT to ' BCT_directory]);
    
    %make parameter file
    fid = fopen(fullfile(outputDir, 'analysis_parameters.txt'),'w');
    fprintf(fid,'List of paramters\n\n');
    fprintf(fid,'\tFrames per second = %f\n', real_FPS);
    fprintf(fid,'\tFC_method = %s\n', params.FC_method);
    fprintf(fid,'\tFC_inactive = %s\n', params.FC_inactive);
    fprintf(fid,'\tFC_islands = %s\n', FC_islands);
    fprintf(fid,'\tCa event threshold = %.2f\n', params.event_thresh);
    fprintf(fid,'\tNo weight threshold applied to FC matrix\n');
    fclose(fid);
    
    %write to separate files for all the output vectors
    output_by_cell = horzcat((1:n)',output.degree_distribution, output.strength_distribution,...
        output.eccentricities, output.clustering_coeff_distribution,...
        output.betweenness_distribution);
    header = ['ROI Number,Degree,Strength,Eccentricity,Clustering Coefficient,'...
        'Betweenness Centrality (normalized)'];
    fid = fopen(fullfile(BCT_directory,'per_cell_measures.csv'),'w');
    fprintf(fid,'%s\n',header);
    fclose(fid);
    dlmwrite(fullfile(BCT_directory,'per_cell_measures.csv'), output_by_cell, '-append');
    
    %rich club curve has to be in its own file bc its indexed differently
    csvwrite(fullfile(BCT_directory, 'rich_club_curve.csv'), output.rich_club_curve);
    
    %write the rest of the output struct to 'network_summary.csv'
    output2 = output;
    output2.degree_distribution = 'written to separate file';
    output2.strength_distribution = 'written to separate file';
    output2.clustering_coeff_distribution = 'written to separate file';
    output2.betweenness_distribution = 'written to separate file';
    output2.rich_club_curve = 'written to separate file';
    output2.eccentricities = 'written to separate file';
    
    temp_table = struct2table(output2);
    writetable(temp_table,fullfile(BCT_directory, 'network_summary.csv'));
    save(fullfile(BCT_directory, 'BCT_output.mat'), 'output');
    
    %save events & cell summary stats
    cell_directory = fullfile(outputDir, 'Events_Cell');
    if ~exist(cell_directory, 'dir')
        mkdir(cell_directory);
    end
    cprintf('*blue','%s\n', ['Saving events & cell metrics to ' cell_directory]);
    events = processed_analysis.dat(:,2);
    single_cell_summary = {
        'Mean Events' mean(events);
        'Mean Events per Second' mean(events)/(processed_analysis.Frames/real_FPS);
        'Mean Events (Active Cells)' mean(events(events > 0));
        'Mean Events per Second (Active Cells)' mean(events(events > 0))/(processed_analysis.Frames/real_FPS);
        'Standard Deviation Events' std(events);
        '0 Bin' numel(events(events == 0));
        '1-3 Bin' numel(events(events > 0 & events < 4));
        '4-6 Bin' numel(events(events > 3 & events < 7));
        '7-9 Bin' numel(events(events > 6 & events < 10));
        '10+ Bin' numel(events(events > 9));
        'Max Events' max(events);
        'Min Events' min(events);
        'Number of Cells' numel(events);
        'Number of Active Cells' numel(events(events > 0));
        'Fraction of Cells That Are Active' numel(events(events > 0))/numel(events);
        };
    
    dlmwrite(fullfile(cell_directory, 'events_cell_summary_stats.csv'),[]);  %clear the file
    fid = fopen(fullfile(cell_directory, 'events_cell_summary_stats.csv'),'wt');
    if fid>0
        for k=1:size(single_cell_summary,1)
            fprintf(fid,'%s,%f\n',single_cell_summary{k,:});
        end
        fclose(fid);
    end
    
    %save cell by cell event data
    dlmwrite(fullfile(cell_directory, 'events_by_cell.csv'),[]); %clear the file
    header = {'ROI Number' 'Total Events' 'Events/Sec'};
    headerJoined = strjoin(header, ',');
    fid = fopen(fullfile(cell_directory, 'events_by_cell.csv'),'w');
    fprintf(fid,'%s\n',headerJoined);
    fclose(fid);
    data = horzcat((1:length(events))', events, events./(processed_analysis.Frames/real_FPS));
    dlmwrite(fullfile(cell_directory, 'events_by_cell.csv'),data,'-append');
    
    %make events per cell histogram
    f = figure();
    histogram(events,max(events));
    xlabel('Events per cell');
    ylabel('Number of cells');
    title('Event Histogram');
    saveas(f,fullfile(cell_directory, 'histogram.fig'));
    
    %computation for Ca event and interspike interval raster plot
    spikes = processed_analysis.Spikes_cell;
    num_active_cells = sum(events > 0);
    
    raster_x = zeros(sum(events),1);
    raster_y = zeros(sum(events),1);
    interspike_x = zeros(sum(events) - num_active_cells, 1);
    interspike_y = zeros(sum(events) - num_active_cells, 1); 
    p = 1;
    interspike_p = 1;
    
    for w = 1:length(spikes)
        for j = 1:length(spikes{w})
            raster_x(p) = spikes{w}(j);
            raster_y(p) = w;
            p = p + 1;
            
            if j > 1
                interspike_x(interspike_p) = spikes{w}(j) - spikes{w}(j-1);
                interspike_y(interspike_p) = w;
                interspike_p = interspike_p + 1;
            end   
        end
    end
    
    %plot Ca event raster plot
    raster_x = raster_x ./ real_FPS;
    f = figure();
    scatter(raster_x,raster_y,15);
    xlabel('Time (s)');
    ylabel('Cell number');
    title('Ca Event Scatter plot');
    saveas(f,fullfile(cell_directory, 'CaEvent_raster.fig'));
    
    %Make dF/F raster plot
    f = figure();
    time = ([1:processed_analysis.Frames] - 1) / 5;
    cells = 1:processed_analysis.N;    
    surf(time,cells,processed_analysis.dF_cell);
    xlim([0 time(end)]);
    ylim([1 cells(end)]);   
    view(2);
    colorbar;
    xlabel('Time (s)');
    ylabel('Cell number');
    title('Ca Trace (dF/F) Raster Plot');
    saveas(f,fullfile(cell_directory, 'CaTrace_raster.fig'));
    
    %make interspike interval raster plot
    interspike_x = interspike_x ./real_FPS;
    f = figure();
    scatter(interspike_x,interspike_y,15);
    xlabel('Time (s)');
    ylabel('Cell number');
    title('Interspike Interval Scatter plot');
    saveas(f,fullfile(cell_directory, 'InterspikeInterval_raster.fig'));
    
    %sort FC matrix according to modular form
    f = figure();
    fc = processed_analysis.FC.CC.C;
    fc = [fc Cj];
    fc = sortrows(fc, length(fc));
    fc = fc(:,1:end-1);
    fc = fc';
    fc = [fc Cj];
    fc = sortrows(fc, length(fc));
    fc = fc(:,1:end-1);
    fc = fc';
    
    %plot sorted fc matrix
    surf(fc)
    view(2);
    colorbar;
    xlabel('Neurons');
    ylabel('Neurons');
    title('Functional Connectivity Matrix');
    saveas(f,fullfile(cell_directory, 'FC_matrix.fig'));

   
    %make table of all Ca Events
    dlmwrite(fullfile(cell_directory, 'all_events.csv'), []);  %clear the file
    fid = fopen(fullfile(cell_directory, 'all_events.csv'),'w');
    fprintf(fid,'ROI, Event timings (s)\n');
    fclose(fid);
    for x = 1:length(spikes)
        labelled_spikes = [x spikes{x}./real_FPS];
        dlmwrite(fullfile(cell_directory, 'all_events.csv'), labelled_spikes, '-append');
    end
    
end

disp(['Completed analysis on:' strjoin(selectedFolders, '\n')]);

end
