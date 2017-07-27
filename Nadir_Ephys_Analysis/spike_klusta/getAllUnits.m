function units = getAllUnits(spikes, stamps, clusters, features, waves, groups, groupnumber,samplerate)

% extracts cluster information for all clusters in a group and puts into a struct

% samplerate = 32000; % samples / second

myClusters = groups(find(groups(:,2)==groupnumber),1);

for i=1:length(myClusters)
    k=find(clusters==myClusters(i));
    mywaves = extractClusterWaves(waves, clusters, myClusters(i));
    [amp maxchan] = maxChannel(mywaves);
    [rpv id] = spikeQuality(myClusters(i),spikes,clusters,features,3,samplerate,0,'');
    myUnits(i)=struct('group',groupnumber,'cluster',myClusters(i),'spikes',spikes(k),'stamps',stamps(k),'waves',mywaves,'maxChannel',maxchan,'RPV',rpv(1,2),'IsoDist',id(1,2));
end

if isempty(myClusters)
    myUnits=[];
end

units = myUnits;

end

