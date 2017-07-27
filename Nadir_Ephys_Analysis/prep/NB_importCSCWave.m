% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 5, 2016

% This function reshapes the CSC Samples 2D Array into a 1D wave format
function wave = NB_importCSCWave (samples)

% Returns # of rows and columns of the samples matrix
% Rows = # of Samples taken per TimeStamp
% Cols = # of TimeStamps
[rows,cols] = size(samples);

% Elements taken columnwise from the samples file 
% Information reshaped into a continuous 1D array that should have wave
% properties
wave = reshape (samples, 1, rows*cols);

end

