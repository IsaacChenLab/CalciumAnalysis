% plot (tvector (1:10), CSC12_Wave (1:10))

ms = 550;
offsetms = -50;
channels = 32;
numstims = 50;

v = CSC32_Wave;

interp_index = find(v > 30000);

%TTLon = find (TTLs==1);

%stim1 = EV_Timestamps(TTLon);

%interp_index = interp1(tvector,1:length(tvector),stim1,'nearest');

offset = (offsetms/1000)*sampleRate;

lengthWave = (ms/1000)*sampleRate;

% plot CSC2_Wave (startstim:(startstim+3200))
% plot CSC2_Wave(startstim:(startstim+3200))
% plot CSC2_Wave(startstim:3968118)
% plot CSC2_Wave
% plot (CSC2_Wave(startstim:(startstim+3200)))

CSD_DataBlock = zeros (channels, lengthWave, numstims);

for i = 1:numstims;
    
    startstim = interp_index (i) + offset;
    endstim = startstim + lengthWave - 1;

% plot (CSC2_Wave(startstim:(startstim+3200)))

CSD_DataBlock (1,:,i) = CSC1_Wave(startstim:endstim);
CSD_DataBlock (2,:,i) = CSC2_Wave(startstim:endstim);
CSD_DataBlock (3,:,i) = CSC3_Wave(startstim:endstim);
CSD_DataBlock (4,:,i) = CSC4_Wave(startstim:endstim);
CSD_DataBlock (5,:,i) = CSC5_Wave(startstim:endstim);
CSD_DataBlock (6,:,i) = CSC6_Wave(startstim:endstim);
CSD_DataBlock (7,:,i) = CSC7_Wave(startstim:endstim);
CSD_DataBlock (8,:,i) = CSC8_Wave(startstim:endstim);
CSD_DataBlock (9,:,i) = CSC9_Wave(startstim:endstim);
CSD_DataBlock (10,:,i) = CSC10_Wave(startstim:endstim);
CSD_DataBlock (11,:,i) = CSC11_Wave(startstim:endstim);
CSD_DataBlock (12,:,i) = CSC12_Wave(startstim:endstim);
CSD_DataBlock (13,:,i) = CSC13_Wave(startstim:endstim);
CSD_DataBlock (14,:,i) = CSC14_Wave(startstim:endstim);
CSD_DataBlock (15,:,i) = CSC15_Wave(startstim:endstim);
CSD_DataBlock (16,:,i) = CSC16_Wave(startstim:endstim);
CSD_DataBlock (17,:,i) = CSC17_Wave(startstim:endstim);
CSD_DataBlock (18,:,i) = CSC18_Wave(startstim:endstim);
CSD_DataBlock (19,:,i) = CSC19_Wave(startstim:endstim);
CSD_DataBlock (20,:,i) = CSC20_Wave(startstim:endstim);
CSD_DataBlock (21,:,i) = CSC21_Wave(startstim:endstim);
CSD_DataBlock (22,:,i) = CSC22_Wave(startstim:endstim);
CSD_DataBlock (23,:,i) = CSC23_Wave(startstim:endstim);
CSD_DataBlock (24,:,i) = CSC24_Wave(startstim:endstim);
CSD_DataBlock (25,:,i) = CSC25_Wave(startstim:endstim);
CSD_DataBlock (26,:,i) = CSC26_Wave(startstim:endstim);
CSD_DataBlock (27,:,i) = CSC27_Wave(startstim:endstim);
CSD_DataBlock (28,:,i) = CSC28_Wave(startstim:endstim);
CSD_DataBlock (29,:,i) = CSC29_Wave(startstim:endstim);
CSD_DataBlock (30,:,i) = CSC30_Wave(startstim:endstim);
CSD_DataBlock (31,:,i) = CSC31_Wave(startstim:endstim);
CSD_DataBlock (32,:,i) = CSC32_Wave(startstim:endstim);

end
% save CSD_Data

dt = (1/sampleRate)*1000

CSD_Data = mean (CSD_DataBlock,3);

save CSD_load.mat CSD_Data;

% cd ../../../../Dropbox/NeuroSurg/Code/CSDplotter-0.1.1/
% CSDplotter
% plot (CSC2_Wave(startstim:(startstim+3200)))