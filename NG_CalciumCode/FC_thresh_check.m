%counts the number of FC edges that exceed threshold for each threshold

load('processed_analysis.mat');
FC = processed_analysis.FC.CC.C;

thresholds = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99];
l = length(thresholds);
counts = zeros(l,1);

for i = 1:l
    bool_mat = FC > thresholds(i);
    counts(i) = sum(bool_mat(:))/2;
end

counts

 
    