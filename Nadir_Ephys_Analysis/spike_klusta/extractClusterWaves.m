function waves = extractClusterWaves(allwaves, allclusters, clusternumber)

spikes = find(allclusters==clusternumber);

for i=1:length(spikes)
    temp(:,:,i) = allwaves(:,:,spikes(i));
end

waves = temp;
end
