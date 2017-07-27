% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 10, 2016

% Loads TT, EV, and CSC data into workspace
NB_loadall_NLX;

% Converts data to best of fixed or floating point format with 15 digits for double and 7 digits for single.
format long g;

% Turns 2D CSC_Samples data and CSC_Timestamp into a 1D vector format
NB_prepdataCSC_32chan;

% Averages CSC Wave data into a CSD_DataBlock
NB_prep_CSDdata_IO_oneTrial;

% Plots 32 channels
NB_plot32_IO;

% Saves figures
%axis([-50 500 (-10*10^4) 1*10^4]);
%NB_saveFigures;