% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

% Create timestamp arrays for different cells within tetrode file.

% ** this script can be improved by aggregating cell timestamp data into a
% struct or multidimensional array


% Initialize cell arrays
Cell_1_0=[];
Cell_1_1=[];
Cell_1_2=[];
Cell_1_3=[];
Cell_1_4=[];
Cell_1_5=[];
Cell_1_6=[];
Cell_1_7=[];
Cell_1_8=[];
Cell_1_9=[];
Cell_1_10=[];

% Create cell with spike sample information for each cluster
for spike=1:numSpikes
    if TT1_CellNumbers(spike)==0
        Cell_1_0=[Cell_1_0 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==1
        Cell_1_1=[Cell_1_1 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==2
        Cell_1_2=[Cell_1_2 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==3
        Cell_1_3=[Cell_1_3 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==4
        Cell_1_4=[Cell_1_4 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==5
        Cell_1_5=[Cell_1_5 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==6
        Cell_1_6=[Cell_1_6 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==7
        Cell_1_7=[Cell_1_7 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==8
        Cell_1_8=[Cell_1_8 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==9
        Cell_1_9=[Cell_1_9 TT1_Timestamps(spike)];
    elseif TT1_CellNumbers(spike)==10
        Cell_1_10=[Cell_1_10 TT1_Timestamps(spike)];
    end
end

% Changes timestamps to units of seconds.
% Cell_1_0=Cell_1_0 / 1000000;
% Cell_1_1=Cell_1_1 / 1000000;
% Cell_1_2=Cell_1_2 / 1000000;
% Cell_1_3=Cell_1_3 / 1000000;
% Cell_1_4=Cell_1_4 / 1000000;
% Cell_1_5=Cell_1_5 / 1000000;
% Cell_1_6=Cell_1_6 / 1000000;
% Cell_1_7=Cell_1_7 / 1000000;
% Cell_1_8=Cell_1_8 / 1000000;
% Cell_1_9=Cell_1_9 / 1000000;
% Cell_1_10=Cell_1_10 / 1000000;


% TT1
%  Attempt to use for loop to create cells
% for cluster = 0:numClusters-1
%     eval(sprintf('Cell_%d = []', cluster));
% end
% 
% for spike = 1:numSpikes
%     for cluster = 0:numClusters-1
%         if TT1_CellNumbers(spike) == cluster
%             Cell_0 = [Cell_0 TT1_Timestamps(spike)]
% %             eval(sprintf('Cell_%d', cluster)) = [eval(sprintf('Cell_%d', cluster)) TT1_Timestamps(spike)];        
%         end
%     end
% end
