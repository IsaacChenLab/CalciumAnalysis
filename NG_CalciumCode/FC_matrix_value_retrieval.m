load('processed_analysis.mat');
w_FC = processed_analysis.FC.CC.C;
uw_FC = processed_analysis.FC.CC.A;

ROIs = [1,4,6,26,27,48,83,94,109,169];
L = length(ROIs);

w_values = zeros(L);
uw_values = zeros(L);

for i = 1:L
    for j = 1:L
        
        idx_i = ROIs(i);
        idx_j = ROIs(j);
        
        if(i ~= j)
            w_values(i,j) = w_FC(idx_i,idx_j);
            uw_values(i,j) = uw_FC(idx_i,idx_j);
        end
        
    end
end
            
