% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

% Plots Peri-Stim Time Histograms for each cell

% Total period of time to be plotted
ms = 2000;

% Desired timebin in miliseconds
timebinLength = 100 * 10^3; % multiplying by 10^3 converts miliseconds to microseconds

% Converts plot length into a microsecond timestamp format (x-axis)
lengthPlot = (ms/1000)*sampleRate;

% Initialize Peri-Stim Raster Plot data matrix
% Row  1 = stim timestamps
% Rows 2+ = spike time stamps
PSRPdata = zeros(length(Cell_1_1)+1,length(stimOnTimestamps));
PSRPdata(1,:) = stimOnTimestamps;

% Time between each stim
interStimInterval = stimOnTimestamps(2)-stimOnTimestamps(1);

% Create Peri-Stim Time Histogram

% Divide interStimInterval into timebins
% Find out if spikes occur within a timebin
% Plot timebins with # of spikes

% Calculate number of timebins and initialize timebinVector
numTimebins = interStimInterval/timebinLength;
timebinVector = [zeros(1,20)];

% Populate timebinVector which contains each timebin
for binCounter = 1:numTimebins
    
    timebinVector(binCounter) = timebinLength * binCounter;
    
end

% Initialize PSTHdata matrix
PSTHdata = zeros(length(Cell_1_1)+1,numTimebins);
PSTHdata(1,:) = timebinVector;

% For each stim, find out if spikes occur within a timebin
% For each stim
for thisStim = 1:numStims
   
    % Look at each timebin in the interStimInterval
    for binCounter = 1:numTimebins
       
        % Look at each spike
        for thisSpike = 1:length(Cell_1_1)
            
            % If the spike falls within the timebin
            timebinStart = stimOnTimestamps(thisStim) + timebinLength*(binCounter-1);
            timebinEnd   = stimOnTimestamps(thisStim) + timebinLength*binCounter;
            if(Cell_1_1(thisSpike) >= timebinStart && Cell_1_1(thisSpike) < timebinEnd)
                tempStamp = Cell_1_1(thisSpike);
                
                % Add the cell's timestamp to the first zero in the timebin
                % column
                for columnValue = 1:length(PSTHdata(:,1))
                
                    if PSTHdata(columnValue,binCounter) == 0
                        PSTHdata(columnValue,binCounter) = tempStamp; 
                        break
                    end
                    
                end
                
            end
        
        end
    end
    
end

figure
hold on

% Plot the PSTH
% x = timebins
x = 0:interStimInterval/numTimebins:interStimInterval-1;
% y = number of spikes in timebin
for thisBin = 1:numTimebins    
    % Find number of spikes (subtract 1 because PSTHdata(1,:) is timebin markers
    y(thisBin) = nnz(PSTHdata(:,thisBin))-1;
end
% Plot
psth = bar(x,y,'histc');
set(psth,'facecolor','k','edgecolor','k')
ax = gca;
ax.YTick = 0:1:max(y);
xlabel('Time (seconds)')
ylabel('Number of spikes')

% Find out if spikes occur within an interStimInterval
% for each stim
for thisStim = 1:numStims
    
    % look at each cell's spikes
    for thisSpike = 1:length(Cell_1_1)
        
        % Calculate if the spike timestamp is within a stim window
        % for the last stim, use interStimInterval to calculate stim window
        if thisStim == numStims;
            if (Cell_1_1(thisSpike) >= stimOnTimestamps(thisStim) && Cell_1_1(thisSpike) < (stimOnTimestamps(thisStim)+interStimInterval) == 1)
                tempStamp = Cell_1_1(thisSpike);
                
                % Find the first zero in the stim column and add the spike to PSRPdata
                for i = 1:length(PSRPdata(:,1))
                   if PSRPdata(i,thisStim) == 0
                       PSRPdata(i,thisStim) = tempStamp;
                       break
                   end   
                end
                
            end
            
        % for all other stims, use stim(i+1) - stim(i) to calculate window
        elseif (Cell_1_1(thisSpike) >= stimOnTimestamps(thisStim) & Cell_1_1(thisSpike) < stimOnTimestamps(thisStim+1)) == 1;
           tempStamp = Cell_1_1(thisSpike); 
           for i = 1:length(PSRPdata(:,1))
               if PSRPdata(i,thisStim) == 0
                   PSRPdata(i,thisStim) = tempStamp;
                   break
               end
           end
        end
        
    end
    
end

% Now, PSRPdata in row 1 contains stimOnTimestamps
% and PSRPdata in row 2 and below contains spike timestamps

figure('units','normalized','outerposition',[0 0 1 1]);
hold on

% Plot the PSRP
for stim = 1:numStims   
    % y value = trial number
    y = stim;
    
    % look at PSRPdata for each stim and see if there are spikes within the
    % window
    for spike = 2:length(Cell_1_1)-1
        if PSRPdata(spike, stim) ~= 0
            x = (PSRPdata(spike,stim) - PSRPdata(1,stim))/1000000; % divide by 1 million = microseconds to seconds
            plot(x,y,'k*')
            xlabel('Time')
            ylabel('Trials')
        end
    end
    
end

axis([0 interStimInterval/1000000 0 numStims])



% Create figure with both graphs

figure('units','normalized','outerposition',[0 0 1 1]);
hold on


subplot(2,1,1)
hold on
% Plot the PSRP
for stim = 1:numStims   
    % y value = trial number
    y = stim;
    
    % look at PSRPdata for each stim and see if there are spikes within the
    % window
    for spike = 2:length(Cell_1_1)-1
        if PSRPdata(spike, stim) ~= 0
            x = (PSRPdata(spike,stim) - PSRPdata(1,stim))/1000000; % divide by 1 million = microseconds to seconds
            plot(x,y,'k*')
%            xlabel('Time')
            ylabel('Trials')
        end
    end
    
end

axis([0 interStimInterval/1000000 0 numStims])
title(['n = ' num2str(nnz(PSTHdata(2:end,:)))])

subplot(2,1,2)
hold on
x = (0:interStimInterval/numTimebins:interStimInterval-1)/1000000;
% y = number of spikes in timebin
for thisBin = 1:numTimebins    
    % Find number of spikes (subtract 1 because PSTHdata(1,:) is timebin markers
    y(thisBin) = nnz(PSTHdata(:,thisBin))-1;
end
% Plot
psth = bar(x,y,'histc');
set(psth,'facecolor','k','edgecolor','k')
ax = gca;
ax.YTick = 0:1:max(y);
xlabel('Time (seconds)')
ylabel('Number of spikes')

