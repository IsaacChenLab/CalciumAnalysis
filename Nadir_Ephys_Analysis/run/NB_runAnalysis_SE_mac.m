% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% June 22, 2016

% Loads TT (tetrode), EV (event), and CSC (continuously sampled channel) data into workspace
NB_loadall_NLX_v3;

% Converts data to best of fixed or floating point format with 15 digits for double and 7 digits for single
% Basically, cuts data off at the 15th decimal place
format long g;

% Turns 2D CSC_Samples data and CSC_Timestamp into a 1D vector format
NB_prepdataCSC_32chan;

% Adds all CSC data into CSD_DataBlock
% Averages CSC Wave data into a CSD_DataBlock
NB_prep_CSDdata_Stim;

% Plots 32 channels
NB_plot32;

% Check responses
%NB_plot32_picker;

% Plot positive and negative responses
%NB_plot32_SE_rateHistogram;

% Adjust axes
%NB_plot_adjust_electrical_subplots


% Saves figures
%NB_saveFigures;