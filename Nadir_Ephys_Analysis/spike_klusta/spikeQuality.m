function [RPVs IDs] = spikeQuality(myClusters, spikes, clusters, features, RPVthresh, samplerate, writedata, filebase)

% Calculates Isolation Distance and Refractory Period Violations and writes
% them to filebase.kwik if writedata = 1

% myClusters - vector size M of cluster numbers of interest as assigned by
% klusta
% spikes - all spike times in dataset
% clusters - cluster assignments of all spikes in dataset
% features - N x D array where N is number of spikes in dataset, D is
%   number of features (3 x number of channels)
% RPVthresh, samplerate - see function RPV
% writedata - 1 to write calculations to kwik file

% Output:
% RPV - M x 3 matrix: first column is cluster number, second column is
% number of refractory period violations, third column is total spikes
%
% ID - M x 2 matrix: first colum is cluster number, second column is
% Isolation Distance

i=1;
RPVtemp = zeros(length(myClusters),3);
IDtemp = zeros(length(myClusters),2);

for j=myClusters'
    RPVtemp(i,1)=j;
    [RPVtemp(i,2) RPVtemp(i,3)] = RPV(spikes(find(clusters==j)),RPVthresh,samplerate);
    IDtemp(i,1)=j;
    IDtemp(i,2) = IsolationDistance_pfk(features,find(clusters==j));
    
    if writedata == 1
        h5writeatt([filebase '.kwik'],['/channel_groups/0/clusters/main/' num2str(j) '/quality_measures'],'refractory_violation',RPVtemp(i,2));
        h5writeatt([filebase '.kwik'],['/channel_groups/0/clusters/main/' num2str(j) '/quality_measures'],'isolation_distance',IDtemp(i,2));
    end
   
    i=i+1;
end



RPVs = RPVtemp;
IDs = IDtemp;
end
