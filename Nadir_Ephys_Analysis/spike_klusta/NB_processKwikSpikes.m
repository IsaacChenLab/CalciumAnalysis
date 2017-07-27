% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

% Setup tvector varible
LoadCSCs_pfk 
tvectors = struct('tvector', tvector);

% Setup datablock variable
myid = fopen('kwikdata.dat', 'r')
rawdatablock = fread(myid, [1, inf], 'int16');
fclose(myid);
datablock = CSCBandPass(rawdatablock, [600 6000], samplerate, 3);

clusterData = processKwikSpikes('kwikdata', 1, 2, tvectors, datablock, 32000);

save clusterData clusterData;

% function mydv = processKwikSpikes(filebase,numdvs,group,tvectors, datablock,samplerate)
% 
% These are the things you need to run processKwikSpikes:
% 
% 1.	the “.kwik” and “.kwx” files created by klusta. These should have the same filebase, i.e. “myclusters.kwik” and “mycluster.kwx”. They need to be in your working directory. Then:
%     a.	filebase = “myclusters”
% 2.	numdvs = the number of different dvs you are going to process
% 3.	group = the number of the cluster group you wish to process (“good” = 2 unless you have changed it)
% 4.	tvectors – this is a little bit of a confusing variable. The is a struct that contains N different tvector variables, where N is the number of dvs you are processing. I am thinking of getting rid of the multiple-dv capability since it is just confusing for most processing, which is just one dv. Anyway, if you are just processing one dv then you need to get the tvector into a struct. Here is an example set of steps:
%     a.	go to the directory that has the CSCs for this recording.
%     b.	Load 1 CSC: LoadCSCs_pfk – this should create a variable called “tvector”
%     c.	Create a new variable of type struct that has your tvector in it:
%         i.	tvectors = struct(‘tvector’, tvector);
% 5.	datablock – this is a matrix of the raw data you fed into klusta. It is contained in your “.dat” file. To create this you need to know how many channels are in your data. The first step is to load this data from the dat file into a variable. The second step is to filter that data for higher frequencies if you would like to have meaningful amplitudes—if you don’t want this then you can use the raw data. These are the steps:
%     a.	myid = fopen(‘myrawdata.dat’, ‘r’);   % this opens the file for reading
%     b.	rawdatablock = fread(myid, [N, inf], ‘int16’);   % this is asking matlab to read in the data into a matrix that has N rows (N = number of channels) and however many columns are needed to read in all data. ‘int16’ is referring to the datatype of each sample (a 16 bit signed integer). This step could take a while.
%     c.	fclose(myid);   % important to close the file when finished reading from it
%     d.	datablock = CSCBandPass(rawdatablock, [600 6000], samplerate, 3);   % this is butterworth filtering the data between 600 and 6000 Hz. Samplerate is usually 32000, or 30000 for wireless data. This step could take quite a long time.
% 6.	samplerate – whatever the sample rate is for your datablock
% 
% So, in this example, you could then run the command:
% 
% processKwikSpikes(‘myclusters’, 1, 2, tvectors, datablock, 32000);
