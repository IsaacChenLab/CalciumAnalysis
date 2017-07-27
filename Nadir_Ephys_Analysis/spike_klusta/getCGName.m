function name = getCGName(clustergroup,filename)

% reads .kwik file and returns the name associated with the cluster group
% number

name = h5readatt(filename,['/channel_groups/0/cluster_groups/main/' num2str(clustergroup)],'name');

end

    