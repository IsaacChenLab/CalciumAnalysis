function [ChenNetworkBatch] = ChenNetworkBatch(real_FPS, only_active_neurons)
%Computes functional network analysis metrics
    % error if no FPS
    if exist('real_FPS', 'var') == 0
        error('Frames per second is not set.');
    end
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
        % Run FluoroSNNAP Processing
        processed_analysis = PostProcess_test(selected_folder);

    %BRAIN CONNECTIVITY TOOLBOX 
        % CROSS CORRELATION OUTPUT
                % Explanation: Matrix of the largest pair-wise correlation (r values) between neurons' firing ...
                % ... across discrete units of lag time.
                    % - This is the matrix output from the FC cross correlation from FluoroSNNAP
                    % - BCT analysis requires an FC matrix of some sort 
                    % - This cross-corr method randomly generates fluorescent traces based on spike times ...
                    % ... then correlates these traces between 2 neurons across discrete units of lag time.
                    % - I don't understand why FluoroSNNAP does not use the raw fluorescence traces.
    disp('Running BCT network analysis');
    %Get FC Matrix to Run BCT
    load([selected_folder slash 'processed_analysis.mat']);
    fc_matrix = processed_analysis.FC.CC.C;
    % Check only active neurons
    if exist('only_active_neurons', 'var') == 0
        only_active_neurons = 0;
    end
    % only active neurons 
    if only_active_neurons == 1
        fc_matrix( all(~fc_matrix,2), : ) = [];
        fc_matrix( :, all(~fc_matrix,1) ) = [];
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
            output.degree_distribution = degrees_und(fc_matrix);
            output.degree_avg = mean(output.degree_distribution); 
            output.degree_avg_normalized = output.degree_avg/size(fc_matrix,1);  
        %CLUSTERING COEFFICIENT
            % Explanation: Average intensity of all triangles associated with each node. (fxnal integration measure)
                % Computed according to Equation 9 defined in Onnela et al. 2005
                % For each cell, it multiplies weights of connections to two other cells, then sums these per cell, ...
                % ... then divides by k(k-1), where k = the degree (defined above)
                % dividing by k(k-1) effectively scales the sum of intensity into an average per cell
            [C_pos,C_neg,output.average_clustering_coeff,Ctot_neg] = clustering_coef_wu_sign(fc_matrix,1);
            output.clustering_coeff_distribution = clustering_coef_wu(fc_matrix);
        %MODULARITY
            [Ci output.modularity] = modularity_und(fc_matrix);
        %CHARACTERISTIC PATH
            % Explanation: Average path between all nodes in network (fxnal segregation measure)
                % A path between two nodes is the inverse of the connectivity weight (1/weight).
                % A shorter path means two nodes are better connected.
                % Therefore, we must first convert the FC matrix into a 'distance' matrix of paths.
                % All the nonzero, nonNaN paths between all nodes are then averaged into one number.
            distance_matrix = weight_conversion(fc_matrix,'lengths');
            D = distance_matrix; % let's get rid of zero columns
            D( all(~D,2), : ) = [];
            D( :, all(~D,1) ) = [];
            [output.chracteristic_path, output.global_efficiency] = charpath(D);
        %GLOBAL EFFICIENCY
            % Explanation: Average inverse of characteristic path for all nodes in the network (fxnal segregation measure)
            %output.global_efficiency = mean(1./D);    
        %BETWEENNESS CENTRALITY
            % Explanation: Fraction of paths in network that contain each node
                % The nodes with greatest centrality are "hub" nodes
                % The algorithm finds the number of all shortest paths between every set of two nodes that ...
                % ...passes through any one node. The greater number of paths through the node, the greater centrality
                % PROBLEM: Our centrality ouptut only shows lots of zero entries...I would expect all values to have ...
                % ...nonzero entries. I will consider reaching out to the algorithm creator.
            B = betweenness_wei(distance_matrix);
            %must normalize
            n = numel(fc_matrix);
            output.betweenness = B./((n-1)*(n-2));
        %save output;
            BCT_directory = strcat(selected_folder, slash, 'BCT');
            mkdir(BCT_directory);
            cprintf('*blue','%s\n', ['Saving BCT to ' BCT_directory]);
            csvwrite(strcat(BCT_directory, slash, 'degree_distribution.csv'), output.degree_distribution);
            csvwrite(strcat(BCT_directory, slash, 'clustering_coeff_distribution.csv'), output.clustering_coeff_distribution);
            csvwrite(strcat(BCT_directory, slash, 'betweenness.csv'), output.betweenness);
            output2 = output;
            output2.degree_distribution = 'double';
            output2.clustering_coeff_distribution = 'double'; 
            output2.betweenness = 'double';
            temp_table = struct2table(output2);
            writetable(temp_table,strcat(BCT_directory, slash, 'network_summary.csv'));
            
         %save events & cell summary stats
            cell_directory = strcat(selected_folder, slash, 'Events_Cell');
            mkdir(cell_directory);
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
            fid = fopen([cell_directory slash 'events_cell_summary_stats.csv'],'wt');
            if fid>0
                for k=1:size(single_cell_summary,1)
                    fprintf(fid,'%s,%f\n',single_cell_summary{k,:});
                end
                fclose(fid);
            end
         %save cell by cell event data
            header = {'ROI Number' 'Total Events' 'Events/Sec'};
            headerJoined = strjoin(header, ',');
            fid = fopen([cell_directory slash 'events_by_cell.csv'],'w');
            fprintf(fid,'%s\n',headerJoined);
            fclose(fid);
            data = horzcat(transpose(1:numel(events)), events, events./(processed_analysis.Frames/real_FPS));
            dlmwrite([cell_directory slash 'events_by_cell.csv'],data,'-append');
    end
          disp(['Completed analysis on:' strjoin(selected_folders, '\n')]);  
end

