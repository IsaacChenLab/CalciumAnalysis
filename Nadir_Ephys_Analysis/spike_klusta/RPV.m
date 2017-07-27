function [violations spiketotal] = RPV(spiketimes, thresh, samplerate)

% returns number of spikes and percentage of spikes that occur within refractory period set by
% thresh (in msec). samplerate is in Hz

tsamples = (samplerate * thresh) / 1000;
count = 0;

for i=1:(size(spiketimes,1) - 1)
    if (spiketimes(i+1) - spiketimes(i)) < tsamples
        count = count + 1;
    end
end

violations = count;
spiketotal = length(spiketimes);
disp 'total spikes';
disp (length(spiketimes));
end
