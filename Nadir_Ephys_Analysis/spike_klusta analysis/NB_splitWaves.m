% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

selectedLine = ginput(2);
line = selectedLine';

goodClusters = zeros(length(clusterData.clusters),1)

for cluster = 1:length(clusterData.clusters);
    if clusterData.clusters(cluster) == 2
        goodClusters(cluster) = 2;
    end
end

figure
hold on
for spike = 1:length(clusterData.spikes);
    if clusterData.clusters(spike) == 0
        plot(clusterData.waves(1,:,spike));
    end
end

for point = 1:32;
    wave(1,point) = point;
end

wave(2,:) = clusterData.waves(1,:,1);

p = InterX(wave, line)


[x0,y0] = intersections(wave(1,:), wave(2,:), line(1,:), line(2,:));