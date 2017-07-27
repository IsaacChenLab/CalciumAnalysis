function [allspikes allstamps allclusters allfeatures clustergroups] = getKwikSpikeData(filebase, tvector)

% gets spike data from a .kwik and .kwx file (output of klusta)
% Reads from filebase.kwik and filebase.kwx
% tvector - list of neuralynx (or other) absolute timestamps for each 
% sample in the dataset - this can be created by running LoadCSCs_pfk for 1 CSC channel

% For N number of spikes:

% allspikes - size N : sample number for each spike (i.e. for a sample rate of 32kHz
% a spike number of 32000 would be 1 second after beginning of dataset)
% allstamps - the absolute time stamp of each spike as referenced by tvector
% allclusters - size N : cluster number as assigned by klusta

% allfeatures - N x D array where N is number of spikes, D is number of
% features (12 for tetrode data: 3 per channel)

% clustergroups - M x 2 array , where M is number of unique clusters; first
% column is the cluster number, second column is the cluster group number
% it belongs to

% use function getCGName to extract text name of cluster group

allspikes = hdf5read([filebase '.kwik'],'/channel_groups/0/spikes/time_samples');
allstamps = tvector(allspikes)';
allclusters = hdf5read([filebase '.kwik'],'/channel_groups/0/spikes/clusters/main');

fmasks = hdf5read([filebase '.kwx'],'/channel_groups/0/features_masks');
allfeatures = (squeeze(fmasks(1,:,:)))';

j=1;
for i=(unique(allclusters))'
    clustergroups(j,1)=i;
    clustergroups(j,2)=h5readatt([filebase '.kwik'],['/channel_groups/0/clusters/main/' num2str(i)],'cluster_group');
    j=j+1;
end

end
