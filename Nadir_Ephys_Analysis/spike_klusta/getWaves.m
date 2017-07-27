function waves = getWaves(spikes, datablock, begin_win, end_win)
% extracts waveforms from data block
% spikes - vector containing the time index (sample number) of detected
% spikes
%
% datablock - N x M matrix of sample values where N is number of channels,
% M is number of samples; will usually need to be filtered in some way to
% eliminate low frequency fluctuations in baseline (i.e. 600-6000 Hz) if
% want to have comparable amplitudes
%
% begin_win - number of samples to extract preceeding spike index
% end_win - number of samples to extract after spike index
%
% waves -  N x P x L matrix  where L is the number of spikes, P is the number of samples per waveform (begin_win + end_win + 1),
% N is the number of channels
% 

    tempwaves = zeros(size(datablock,1),(begin_win+end_win+1),length(spikes));

    for i=1:length(spikes)
        for j=1:size(datablock,1)
            sampleindex = [(spikes(i) - begin_win):(spikes(i) + end_win)];
            if (spikes(i) - begin_win) < 1
                sampleindex = [1:(spikes(i)+end_win)];
            end
            if (spikes(i) + end_win) > size(datablock,2)
                sampleindex = [(spikes(i) - begin_win):size(datablock,2)];
                disp(sampleindex);
            end
            tempwaves(j,1:length(sampleindex),i)=datablock(j,sampleindex);
        end
    end
    
    waves = tempwaves;
end
