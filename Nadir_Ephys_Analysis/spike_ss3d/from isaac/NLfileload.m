
% Load Events.nev file and multiple NTT files.
% x=NLfileload(pathfolder, NTTnumber)
% "Pathfolder" is location of files.
% "NTTnumber" is the number of NTT files to be loaded.

% pathfolder='C:\Users\admin\Desktop\MEAstretch_ephys\091414';
% % NTTnumber=2;
% 
% Eventpath=strcat(pathfolder, '\Events.nev');
% [TimestampEvents, TTLs]=Nlx2MatEV(Eventpath, [1 0 1 0 0], 0, 1, []);
% 
% [TimestampsTT1, CellNumbersTT1]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT1_cells.NTT', [1 0 1 0 0], 0, 1, []);
% [TimestampsTT2, CellNumbersTT2]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT2_cells.NTT', [1 0 1 0 0], 0, 1, []);
% [TimestampsTT3, CellNumbersTT3]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT3_cells.NTT', [1 0 1 0 0], 0, 1, []);
% [TimestampsTT4, CellNumbersTT4]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT4_cells.NTT', [1 0 1 0 0], 0, 1, []);
% [TimestampsTT5, CellNumbersTT5]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT5_cells.NTT', [1 0 1 0 0], 0, 1, []);
% [TimestampsTT6, CellNumbersTT6]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT6_cells.NTT', [1 0 1 0 0], 0, 1, []);
% [TimestampsTT7, CellNumbersTT7]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT7_cells.NTT', [1 0 1 0 0], 0, 1, []);
% [TimestampsTT8, CellNumbersTT8]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT8_cells.NTT', [1 0 1 0 0], 0, 1, []);

%Nadir
% Create arrays for timestamps, cell numbers (groups in SS3D), and samples
[TT1_Timestamps, TT1_CellNumbers, TT1_Samples] = Nlx2MatSpike('TT1_cells.NTT', [1 0 1 0 1], 0, 1, []);

% Create arrays for event timestamps and TTL information
[EV_Timestamps, EV_TTLs] = Nlx2MatEV('Events.nev', [1 0 1 0 0], 0, 1, []);

% Convert 3D TT1 data into 2D data (only use 1 of 4 tetrode channels)
SamplesCSC1 = squeeze(TT1_Samples(:,1,:));

spikeData = struct('TimestampsTT1', TT1_Timestamps, 'CellNumbersTT1', TT1_CellNumbers, 'SamplesTT1', TT1_Samples, 'EV_Timestamps', EV_Timestamps, 'EV_TTLs', EV_TTLs, 'SamplesCSC1', SamplesCSC1);

numSpikes=length(TT1_CellNumbers);
numClusters=length(unique(TT1_CellNumbers));

save spikeData spikeData;

%/Nadir


% load second data set
% pathfolder2='C:\Users\admin\Desktop\MEAstretch_ephys\091414';
% 
% Eventpath2=strcat(pathfolder2, '\Events.nev');
% [TimestampEvents2, TTLs2]=Nlx2MatEV(Eventpath2, [1 0 1 0 0], 0, 1, []);
% 
% [TimestampsTT52, CellNumbersTT52]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT5_cells.NTT', [1 0 1 0 0], 0, 1, []);
% 
% % load third data set
% pathfolder3='C:\Users\admin\Desktop\MEAstretch_ephys\021514';
% 
% Eventpath3=strcat(pathfolder3, '\Events.nev');
% [TimestampEvents3, TTLs3]=Nlx2MatEV(Eventpath3, [1 0 1 0 0], 0, 1, []);
% 
% [TimestampsTT23, CellNumbersTT23]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\021514\TT2_cells.NTT', [1 0 1 0 0], 0, 1, []);
% 

% % load fourth data set
% pathfolder4='C:\Users\admin\Desktop\MEAstretch_ephys\091414';
% 
% Eventpath4=strcat(pathfolder4, '\Events.nev');
% [TimestampEvents4, TTLs4]=Nlx2MatEV(Eventpath4, [1 0 1 0 0], 0, 1, []);
% 
% [TimestampsTT54, CellNumbersTT54]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\091414\TT5_cells.NTT', [1 0 1 0 0], 0, 1, []);
% 
% % load fifth data set
% pathfolder5='C:\Users\admin\Desktop\MEAstretch_ephys\122813';
% 
% Eventpath5=strcat(pathfolder5, '\Events.nev');
% [TimestampEvents5, TTLs5]=Nlx2MatEV(Eventpath5, [1 0 1 0 0], 0, 1, []);
% 
% [TimestampsTT85, CellNumbersTT85]=Nlx2MatSpike('C:\Users\admin\Desktop\MEAstretch_ephys\122813\TT8_cells.NTT', [1 0 1 0 0], 0, 1, []);




% Attempt to automate & generalize to any number of NTT files.

% for a=1:NTTnumber
%     Timestamps.(genvarname(['TimestampsTT', num2str(a)]))=0;
%     ScNumbers.(genvarname(['ScNumbersTT', num2str(a)]))=0;
%     CellNumbers.(genvarname(['CellNumbersTT', num2str(a)]))=0;
%     Features.(genvarname(['FeaturesTT', num2str(a)]))=0;
%     Samples.(genvarname(['SamplesTT', num2str(a)]))=0;
%     
%     NTTpath=strcat(pathfolder, '\TT', num2str(a), '_cells.NTT');
%     
%     [Timestamps, ScNumbers, CellNumbers, Features, Samples]=Nlx2MatSpike(NTTpath, [1 1 1 1 1], 0, 1, []);
% end

% for a=1:NTTnumber
%     Timestamps=strcat('TimestampsTT', num2str(a));
%     ScNumbers=strcat('ScNumbersTT', num2str(a));
%     CellNumbers=strcat('CellNumbersTT', num2str(a));
%     Features=strcat('FeaturesTT', num2str(a));
%     Samples=strcat('SamplesTT', num2str(a));
%     
%     NTTpath=strcat(pathfolder, '\TT', num2str(a), '_cells.NTT');
%     
%     [Timestamps, ScNumbers, CellNumbers, Features, Samples]=Nlx2MatSpike(NTTpath, [1 1 1 1 1], 0, 1, []);
% end

