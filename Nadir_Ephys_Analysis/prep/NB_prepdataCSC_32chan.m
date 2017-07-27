% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 5, 2016

% This file prepares CSC data and Timestamp data into vector formats

% Convert CSC_Sample data into vector format
CSC1_Wave = NB_importCSCWave (CSC1_Samples);
CSC2_Wave = NB_importCSCWave (CSC2_Samples);
CSC3_Wave = NB_importCSCWave (CSC3_Samples);
CSC4_Wave = NB_importCSCWave (CSC4_Samples);
CSC5_Wave = NB_importCSCWave (CSC5_Samples);
CSC6_Wave = NB_importCSCWave (CSC6_Samples);
CSC7_Wave = NB_importCSCWave (CSC7_Samples);
CSC8_Wave = NB_importCSCWave (CSC8_Samples);
CSC9_Wave = NB_importCSCWave (CSC9_Samples);
CSC10_Wave = NB_importCSCWave (CSC10_Samples);
CSC11_Wave = NB_importCSCWave (CSC11_Samples);
CSC12_Wave = NB_importCSCWave (CSC12_Samples);
CSC13_Wave = NB_importCSCWave (CSC13_Samples);
CSC14_Wave = NB_importCSCWave (CSC14_Samples);
CSC15_Wave = NB_importCSCWave (CSC15_Samples);
CSC16_Wave = NB_importCSCWave (CSC16_Samples);
CSC17_Wave = NB_importCSCWave (CSC17_Samples);
CSC18_Wave = NB_importCSCWave (CSC18_Samples);
CSC19_Wave = NB_importCSCWave (CSC19_Samples);
CSC20_Wave = NB_importCSCWave (CSC20_Samples);
CSC21_Wave = NB_importCSCWave (CSC21_Samples);
CSC22_Wave = NB_importCSCWave (CSC22_Samples);
CSC23_Wave = NB_importCSCWave (CSC23_Samples);
CSC24_Wave = NB_importCSCWave (CSC24_Samples);
CSC25_Wave = NB_importCSCWave (CSC25_Samples);
CSC26_Wave = NB_importCSCWave (CSC26_Samples);
CSC27_Wave = NB_importCSCWave (CSC27_Samples);
CSC28_Wave = NB_importCSCWave (CSC28_Samples);
CSC29_Wave = NB_importCSCWave (CSC29_Samples);
CSC30_Wave = NB_importCSCWave (CSC30_Samples);
CSC31_Wave = NB_importCSCWave (CSC31_Samples);
CSC32_Wave = NB_importCSCWave (CSC32_Samples);

% Convert CSC_Timestamp data into vector format
timeVector = NB_timeVector3(CSC1_Timestamp);

% Calculate samplerate
sampleRate = CSC1_SampleFrequency(1);