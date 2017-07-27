% Total period of time to be plotted
ms = 1000;

% Number of miliseconds to include before the stimulation
offsetms = -50;

% Converts offsetms into a microsecond timestamp format (y-axis)
offset = (offsetms/1000)*sampleRate;

% Converts plot length into a microsecond timestamp format (x-axis)
lengthPlot = (ms/1000)*sampleRate;

% Number of channels
channels = 32;

% Find all TTL events when stim is turned on
TTLon = find (EV_TTLs==1);

% Finds all timestamps when stim is turned on based on TTL events
stimOnTimestamps = EV_Timestamps(TTLon);

% Number of stims calculated by number of times stim was turned on
numStims = length(stimOnTimestamps);

spikes = data.stamps';

spikesInWindow = null(length(spikes));

for thisStim = 1:numStims;

    for thisSpike = 1:length(spikes);
        
        if thisStim == numStims;
        elseif (spikes(thisSpike) >= stimOnTimestamps(thisStim) & spikes(thisSpike) < stimOnTimestamps(thisStim+1)) == 1;
            % turn spike stamp into ms value
           tempStamp = spikes(thisSpike) - stimOnTimestamps(thisStim); 
           msStamp = (tempStamp)/1000;
           spikesInWindow(thisSpike) = msStamp;
        end
    end
end

figure;
hold on
for thisSpike = 1:length(spikesInWindow);
    x = spikesInWindow(thisSpike);
    y = 1;
    plot(x,y,'b*');
end


% for rest of time,
% take 2 sec timebins after stim and see what stims look like
timeRemaining = (EV_Timestamps(end) - stimOnTimestamps(numStims))/1000000;

remainingWindowStamps(1) = stimOnTimestamps(end)+2000000;
window = 2;

while remainingWindowStamps(end) < EV_Timestamps(end);
    remainingWindowStamps(window) = remainingWindowStamps(window-1)+2000000;
    window = window+1;
end

for thisWindow = 1:length(remainingWindowStamps);

    for thisSpike = 1:length(spikes);
        
        if thisWindow == length(remainingWindowStamps);
        elseif (spikes(thisSpike) >= remainingWindowStamps(thisWindow) & spikes(thisSpike) < remainingWindowStamps(thisWindow+1)) == 1;
           % turn spike stamp into ms value
           tempStamp = spikes(thisSpike) - remainingWindowStamps(thisWindow); 
           msStamp = (tempStamp)/1000;
           spikesOutOfWindow(thisSpike) = msStamp;
        end
    end
end

for thisSpike = 1:length(spikesOutOfWindow);
    x = spikesOutOfWindow(thisSpike);
    y = 1.1;
    plot(x,y,'r*');
end