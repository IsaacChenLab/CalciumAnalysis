% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

figure;
hold on

for i = 1: length(clusterData.waves); 
    plot(clusterData.waves(1,:,i)); 
end