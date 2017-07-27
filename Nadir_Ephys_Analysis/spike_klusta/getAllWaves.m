function waves = getAllWaves(filename)

% This is now outdated, but can be used for older kwx files that contain
% the actual waveform data.
% Use getWaves instead, which requires a block of raw data

    waves=hdf5read(filename,'/channel_groups/0/waveforms_filtered');
end
