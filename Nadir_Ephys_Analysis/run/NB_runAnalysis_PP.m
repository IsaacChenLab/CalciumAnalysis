% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% June 22, 2016

% Loads TT, EV, and CSC data into workspace
NB_loadall_NLX;

% Converts data to best of fixed or floating point format with 15 digits for double and 7 digits for single.
format long g;

% Turns 2D CSC_Samples data and CSC_Timestamp into a 1D vector format
NB_prepdataCSC_32chan;

% Averages CSC Wave data into a CSD_DataBlock
NB_prep_CSDdata_PP;

% Plots 32 channels
NB_plot32_PP_slopes;

% Saves figures
%axis([-10 100 -4*10^4 1*10^4]);
%NB_saveFigures;