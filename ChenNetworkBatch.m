function ChenNetworkBatch(real_FPS, FC_method)
%Computes functional network analysis metrics
    
    % normally on args are real_fps and only_active_neurons
    % error if no FPS
    
    if exist('real_FPS', 'var') == 0
        error('Frames per second is not set.');
        return;
    end
    
    load('params.mat');
    params.fps = real_FPS;
    fprintf("Ca event threshold = %.2f\n", params.event_thresh);

    if exist('FC_method', 'var') ~= 0
        params.FC_method = FC_method;
    else
         params.FC_method = 'default';
    end
    save('params.mat', '-append');

    
    if strcmpi(params.FC_method, 'raw')
        fprintf("Applying raw trace method for computing FC\n");
    else
        fprintf("Applying FluoroSNNAP default method for computing FC\n");
    end
    fprintf("No threshold for weighted FC matrix being applied\n");
    fprintf("No exclusion of cells without events in FC matrix\n");
    

    
    % handle slashes
    if(ispc)
        slash = '\';
    else
        slash = '/';
    end
    
    % GUI to get folder with .tif file
    selected_folders = uigetfile_n_dir('','Select folders containing tiff stacks');

    % make sure all folders have .tif files
    for i=1:numel(selected_folders)
        selected_folder = selected_folders{i};
        if ~any(size(dir([selected_folder slash '*.tif' ]),1))
            msg = ['There is no .tif file in the selected directory: ' selected_folder];
            error(msg);
            return;
        end
    end
    
    %loop through files & analyze
    disp(['These folders have been selected:' strjoin(selected_folders, '\n')]);  
    for i=1:numel(selected_folders)
        
        selected_folder = selected_folders{i};
        disp(['Batch analysis begun on ' selected_folder]);

        %just check .tif existence again in case someone moved the files
        if ~any(size(dir([selected_folder slash '*.tif' ]),1))       
            msg = 'There is no .tif file in the selected directory.';
            error(msg);
            return;
        end

    % FLUOROSNNAP
        % Run FluoroSNNAP ROI Analysis
        addpath('FluoroSNNAP_code');
        addpath('textprogressbar');
        addpath(['FluoroSNNAP_code' slash 'oopsi-master']);
        AnalyzeData_2(selected_folder, real_FPS);
        fprintf("%s", "check"); 
        % Run FluoroSNNAP Processing
        processed_analysis = PostProcess_test(selected_folder);

        disp('Running BCT network analysis');
        
    %BRAIN CONNECTIVITY TOOLBOX 
        % CROSS CORRELATION OUTPUT
            % Explanation: Matrix of the largest pair-wise correlation (r values) between neurons' firing ...
            % ... across discrete units of lag time.
            % - This is the matrix output from the FC cross correlation from FluoroSNNAP
            % - BCT analysis requires an FC matrix of some sort 
            % - This cross-corr method randomly generates fluorescent traces based on spike times ...
            % ... then correlates these traces between 2 neurons across discrete units of lag time.
            % - I don't understand why FluoroSNNAP does not use the raw fluorescence traces.
        
        %Get FC Matrix to Run BCT
        load([selected_folder slash 'processed_analysis.mat']);
        fc_matrix = processed_analysis.FC.CC.C;
        
        for k=1:2
            switch k 
                case 1
                    disp('BCT network analysis - all neurons');
                    bct_all_vs_active = 'all.';
                case 2
                    disp('BCT network analysis - only active neurons');
                    fc_matrix( all(~fc_matrix,2), : ) = [];
                    fc_matrix( :, all(~fc_matrix,1) ) = [];
                    bct_all_vs_active = 'active.';
            end

        % add BCT to MATLAB path
            addpath('BCT_code');
        %initialize variables
             % we need this for path & centrality calculations
        %metrics of interest
            output = struct();

        %DENSITY
            % Explanation: Density is the fraction of present connections to possible connections.
                % note --> it's all about nonzero FC matrix elements
                % # present connections = # nonzero matrix elements in top triangle of connection matrix 
                % # possible connections = ((matrix size)^2 - (matrix size))/2
                    % note: cannot count self-connections (that's why matrix size is subtracted)
                % vertices => matrix size
                % edges => # nonzero matrix elements in top triangle of connection matrix
            % Execution
            [output.density,output.vertices,output.edges] = density_und(fc_matrix);
        
         %DEGREE DISTRIBUTION
            % Explanation: Degree per cell is the number of connections it has within the network
                % Simply, this is the sum of pair-wise connections (FC matrix elements) between Cell 1 and all ... 
                % ... other cells that are not zero.
                % This would be sensitive to a threshold we could set in generating the FC matrix.
                % Strength is the same as degree but accounts for connection weights
            output.degree_distribution = (degrees_und(fc_matrix))';
            output.degree_avg = mean(output.degree_distribution); 
            output.degree_avg_normalized = output.degree_avg/size(fc_matrix,1);
            output.strength_distribution = (strengths_und(fc_matrix))';
            output.strength_mean = mean(output.strength_distribution);
        
         %CLUSTERING COEFFICIENT
            %Explanation: Clustering coeff of node n is the number of complete triangles there are centered on n
                %weighted by the geometric mean of the edges involved. A complete triangle centered on n is when
                %two of n's neighbors are neighbors with each other. Network average is the average of each node's
                %clustering coeff. The transitivity is the sum of the geometric mean of all triangles in the whole 
                %network divided by the total number of potential triangles (pairs of mutual friends) which is the
                %sum of (degree(n) * degree(n-1)) summed over every node n                
            output.clustering_coeff_distribution = clustering_coef_wu(fc_matrix);
            output.network_avg_clustering_coeff = mean(output.clustering_coeff_distribution);
            output.transitivity = transitivity_wu(fc_matrix);
            %[C_pos,C_neg,output.average_clustering_coeff,Ctot_neg] = clustering_coef_wu_sign(fc_matrix,1);
            
         %MODULARITY
            %Explanation: Modularity of graph is the degree to which that graph can be divided into a set of communities (modules)
                %such that within-community edges are maximized and between-community edges are minimized. Finding the
                % the modularity of the graph requires finding the optimal community structure (ie the one that maximizes
                %modularity) which has been proven to be NP-hard. So this algorithm provides a best estimate of the true
                %modularity. Cj is the optimal community structure (ie the actual modules).
            [Ci, output.modularity] = modularity_und(fc_matrix);

        %LENGTH MATRIX    
            % we must convert the FC matrix into a 'length' matrix of paths.
            % A path between two nodes is the inverse of the connectivity weight (1/weight).
            % A shorter path means two nodes are better connected.
            length_mat = weight_conversion(fc_matrix,'lengths');
            
           
            
            %D = length_matrix;     % let's get rid of zero columns
            %D( all(~D,2), : ) = [];
            %D( :, all(~D,1) ) = [];

            
            
         %DISTANCE MATRIX:
            %Each matrix entry (i,j) is the length of the shortest path between nodes i and j
                %NOTE: Previously the length matrix (where edge weights are inversed so as to represent lengths)
                %was mistaken for the distance matrix. Bug fixed on 7/10/17 by NG
            dist_mat = distance_wei(length_mat);

        %GLOBAL EFFICIENCY
            % Explanation: Average of inverse of distance between all pairs of nodes in the network (fxnal segregation measure)
            %output.global_efficiency = mean(1./D);
        %CHARACTERISTIC PATH
            % Explanation: Average distance between all pairs of nodes in network (fxnal segregation measure)
        %NODE ECCENTRICITY
            % Eccentricity of node n is distance from n to the node that is farthest from n
        %RADIUS
            % the minimum node eccentricity in the network
        %DIAMETER
            % the maximum node eccentricity in the network
            [output.characteristic_path,output.global_efficiency,output.eccentricities,output.radius,output.diameter] = charpath(dist_mat);
            output.mean_eccentricity = mean(output.eccentricities);
            
        %ASSORTATIVITY
            % correlation coefficient for the strengths of nodes that are connected
            output.assortativity = assortativity_wei(fc_matrix, 0); % 0 flag = undirected graph

        %RICH CLUB CURVE
            % for all possible degrees k, the fraction of edges that exist between all nodes of degree k
            % out of the maximum number of such edges that could possibly exist
            % returns a vector of length k where k is the maximum degree of any node in the network
            output.rich_club_curve = (rich_club_wu(fc_matrix))';
            
        %BETWEENNESS CENTRALITY
            % Explanation: Fraction of paths in network that contain each node
                % The nodes with greatest centrality are "hub" nodes
                % The algorithm finds the number of all shortest paths between every set of two nodes that ...
                % ...passes through any one node. The greater number of paths through the node, the greater centrality
            B = betweenness_wei(length_mat);
            %must normalize
            n = length(fc_matrix);
            output.betweenness_distribution = B./((n-1)*(n-2));
            
        %SAVE THE OUTPUT;
            %make the directory
            BCT_directory = strcat(selected_folder, slash, 'BCT');
            if ~exist(BCT_directory, 'dir')
                mkdir(BCT_directory);
            end
            cprintf('*blue','%s\n', ['Saving BCT to ' BCT_directory]);
            
            %write to separate files for all the output vectors
            output_by_cell = horzcat((1:n)',output.degree_distribution, output.strength_distribution, output.eccentricities, output.clustering_coeff_distribution, output.betweenness_distribution);
            header = 'ROI Number,Degree,Strength,Eccentricity,Clustering Coefficient,Betweenness Centrality (normalized)';
            fid = fopen(strcat(BCT_directory,slash,bct_all_vs_active,'per_cell_measures.csv'),'w');
            fprintf(fid,'%s\n',header);
            fclose(fid);
            dlmwrite(strcat(BCT_directory,slash,bct_all_vs_active,'per_cell_measures.csv'), output_by_cell, '-append');

            %rich club curve has to be in its own file bc its indexed differently
            csvwrite(strcat(BCT_directory, slash, bct_all_vs_active, 'rich_club_curve.csv'), output.rich_club_curve);
            
            %write the rest of the output struct to 'network_summary.csv'
            output2 = output;
            output2.degree_distribution = 'written to separate file';
            output2.strength_distribution = 'written to separate file';
            output2.clustering_coeff_distribution = 'written to separate file'; 
            output2.betweenness_distribution = 'written to separate file';
            output2.rich_club_curve = 'written to separate file';
            output2.eccentricities = 'written to separate file';
            
            temp_table = struct2table(output2);
            writetable(temp_table,strcat(BCT_directory, slash, bct_all_vs_active, 'network_summary.csv'));
        end
        
        %save events & cell summary stats
        cell_directory = strcat(selected_folder, slash, 'Events_Cell');
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
    
        dlmwrite([cell_directory slash 'events_cell_summary_stats.csv'],[]);  %clear the file
        fid = fopen([cell_directory slash 'events_cell_summary_stats.csv'],'wt');
        if fid>0
            for k=1:size(single_cell_summary,1)
                fprintf(fid,'%s,%f\n',single_cell_summary{k,:});
            end
            fclose(fid);
        end
        
        %save cell by cell event data
        dlmwrite([cell_directory slash 'events_by_cell.csv'],[]); %clear the file
        header = {'ROI Number' 'Total Events' 'Events/Sec'};
        headerJoined = strjoin(header, ',');
        fid = fopen([cell_directory slash 'events_by_cell.csv'],'w');
        fprintf(fid,'%s\n',headerJoined);
        fclose(fid);
        data = horzcat((1:length(events))', events, events./(processed_analysis.Frames/real_FPS));
        dlmwrite([cell_directory slash 'events_by_cell.csv'],data,'-append');
        
        %make events per cell histogram
        f = figure('Visible', 'off');
        histogram(events,max(events));
        xlabel('Events per cell');
        ylabel('Number of cells');
        title('Event Histogram');
        saveas(f,strcat(cell_directory, slash, 'histogram.jpg'));
        
        %make Ca event raster plot
        spikes = processed_analysis.Spikes_cell;
        x = zeros(sum(events),1);
        y = zeros(sum(events),1);
        p = 1;
        for i = 1:length(spikes)
            for j = 1:length(spikes{i})
                x(p) = spikes{i}(j);
                y(p) = i;
                p = p + 1;
            end
        end
        x = x ./real_FPS;
        f = figure('Visible', 'off');
        scatter(x,y,15,'filled');
        xlabel('Time (s)');
        ylabel('Cell number');
        title('Ca Event Scatter plot');
        saveas(f,strcat(cell_directory, slash, 'raster_plot.jpg'));
        
        %make table of all Ca Events
        dlmwrite(strcat(cell_directory, slash, 'all_events.csv'), []);  %clear the file
        fid = fopen(strcat(cell_directory, slash, 'all_events.csv'),'w');
        fprintf(fid,'ROI, Event timings (s)\n');
        fclose(fid);
        for x = 1:length(spikes)
            labelled_spikes = [x spikes{x}./real_FPS];
            dlmwrite(strcat(cell_directory, slash, 'all_events.csv'), labelled_spikes, '-append');
        end
        
    end
    
  disp(['Completed analysis on:' strjoin(selected_folders, '\n')]);
  
end

