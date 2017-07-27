% plot (tvector (1:10), CSC12_wave (1:10))

offsetgraph = 500;

figure('units','normalized','outerposition',[0 0 1 1]);

hold on

for i = 1:20;
    
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_Data(i,:)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_Data (i,:) - i*offsetgraph;
    
    h = plot(x,y);

hold on;

% CSD_DataBlock (1,:,i) = CSC1_wave(startstim:endstim);
% CSD_DataBlock (2,:,i) = CSC2_wave(startstim:endstim);
% CSD_DataBlock (3,:,i) = CSC3_wave(startstim:endstim);
% CSD_DataBlock (4,:,i) = CSC4_wave(startstim:endstim);
% CSD_DataBlock (5,:,i) = CSC5_wave(startstim:endstim);
% CSD_DataBlock (6,:,i) = CSC6_wave(startstim:endstim);
% CSD_DataBlock (7,:,i) = CSC7_wave(startstim:endstim);
% CSD_DataBlock (8,:,i) = CSC8_wave(startstim:endstim);
% CSD_DataBlock (9,:,i) = CSC9_wave(startstim:endstim);
% CSD_DataBlock (10,:,i) = CSC10_wave(startstim:endstim);
% CSD_DataBlock (11,:,i) = CSC11_wave(startstim:endstim);
% CSD_DataBlock (12,:,i) = CSC12_wave(startstim:endstim);
% CSD_DataBlock (13,:,i) = CSC13_wave(startstim:endstim);
% CSD_DataBlock (14,:,i) = CSC14_wave(startstim:endstim);
% CSD_DataBlock (15,:,i) = CSC15_wave(startstim:endstim);
% CSD_DataBlock (16,:,i) = CSC16_wave(startstim:endstim);
% CSD_DataBlock (17,:,i) = CSC17_wave(startstim:endstim);
% CSD_DataBlock (18,:,i) = CSC18_wave(startstim:endstim);
% CSD_DataBlock (19,:,i) = CSC19_wave(startstim:endstim);
% CSD_DataBlock (20,:,i) = CSC20_wave(startstim:endstim);
% CSD_DataBlock (21,:,i) = CSC21_wave(startstim:endstim);
% CSD_DataBlock (22,:,i) = CSC22_wave(startstim:endstim);
% CSD_DataBlock (23,:,i) = CSC23_wave(startstim:endstim);
% CSD_DataBlock (24,:,i) = CSC24_wave(startstim:endstim);
% CSD_DataBlock (25,:,i) = CSC25_wave(startstim:endstim);
% CSD_DataBlock (26,:,i) = CSC26_wave(startstim:endstim);
% CSD_DataBlock (27,:,i) = CSC27_wave(startstim:endstim);
% CSD_DataBlock (28,:,i) = CSC28_wave(startstim:endstim);
% CSD_DataBlock (29,:,i) = CSC29_wave(startstim:endstim);
% CSD_DataBlock (30,:,i) = CSC30_wave(startstim:endstim);
% CSD_DataBlock (31,:,i) = CSC31_wave(startstim:endstim);
% CSD_DataBlock (32,:,i) = CSC32_wave(startstim:endstim);

end
% save CSD_Data

dt = (1/sampleRate)*1000

% CSD_Data = mean (CSD_DataBlock,3);
% 
% save CSD_load.mat CSD_Data;
% 
% cd ../../../../Dropbox/NeuroSurg/Code/CSDplotter-0.1.1/
% CSDplotter
% plot (CSC2_wave(startstim:(startstim+3200)))