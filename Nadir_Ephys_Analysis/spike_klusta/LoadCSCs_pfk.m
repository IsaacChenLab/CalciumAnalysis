
% load CSCs
% will create a struct with fields: OrigChannel, Timestamp, SampleFrequency, Samples, Wave for
% each selected CSC to load and store them in the array CSC. Thus,
% accessing the Timestamps of CSC14 that was the 5th CSC to be loaded:
% CSC(5).Timestamp, where CSC(5).OrigChannel==14.
% Reformats and preps CSC data into linear arrays called waves (CSC(i).Wave).
% Defines global sample rate and "tvector" based on first CSC loaded.
% CSC.

% If desired, will also load event file, replacing "LoadEVs_v3.m"

clear CSC;      % get rid of any previously loaded data

reply = input('Enter CSC numbers desired to be analyzed, separated by space.\n Leave blank if already defined in array "desiredCSCs".\n >>>','s');
if isempty(reply)
    disp 'Using previously defined CSC list in array "desiredCSCs".'
else
    desiredCSCs = [str2num(reply)];
end
disp('CSCs to use:')
disp (desiredCSCs)
desiredTimestamp = input('Enter full beginning timestamp >>>');
platform = input('Enter 1 if using PC, 0 if mac >>>');
if input('Enter 1 if want to load events >>>')==1
    qualifier=input('Enter events suffix, if any >>>','s');
    
    % Fields for EV loading
    FieldSelection(1) = 1;
    FieldSelection(2) = 1;
    FieldSelection(3) = 1;
    FieldSelection(4) = 0; %Extras
    FieldSelection(5) = 1;

    ExtractHeader = 0;
         
    ExtractMode = 1;
    
    if isempty(desiredTimestamp)
        if isempty(qualifier)
            eventfilename='Events.nev';
        else
            eventfilename=dir(['Events_' qualifier '*.nev']);
        end
    else
        if isempty(qualifier)
            eventfilename=dir(['Events_' num2str(desiredTimestamp) '*.nev']);
            disp 'hello'
        else
            eventfilename=dir(['Events_' qualifier '_' num2str(desiredTimestamp) '*.nev']);
            disp 'hi'
        end
    end
    if platform==1
        [EV_Timestamps, EV_EventIDs, EV_TTLs, EV_EventStrings] = Nlx2MatEV(eventfilename, FieldSelection, ExtractHeader,ExtractMode);
    else
        [EV_Timestamps, EV_EventIDs, EV_TTLs, EV_EventStrings] = Nlx2MatEV_v3(eventfilename.name, FieldSelection, ExtractHeader,ExtractMode);
    end
end
for i = 1:size(desiredCSCs,2)
        if isempty(desiredTimestamp)
            tempfilename=dir(['CSC' num2str(desiredCSCs(i)) '.ncs']);
        else
            tempfilename=dir(['CSC' num2str(desiredCSCs(i)) '_' num2str(desiredTimestamp) '*.ncs']);
        end
        if platform==1
            [tempTimestamp,tempSampleFrequency,tempSamples] = importCSC_v3_PC(tempfilename.name);
        else
            [tempTimestamp,tempSampleFrequency,tempSamples] = importCSC_v3(tempfilename.name);
        end
        CSC(i) = struct('OrigChannel',desiredCSCs(i),'Timestamp',tempTimestamp,'SampleFrequency',tempSampleFrequency,'Samples',tempSamples,'Wave',importCSC_JAW (tempSamples));
        progress = ['CSC',num2str(desiredCSCs(i)),' loaded and transformed to wave.'];
        disp(progress)
end

tvector = NB_timeVector3 (CSC(1).Timestamp);

if CSC(1).SampleFrequency(1)==30303
	samplerate = CSC(1).SampleFrequency(1)./subsample;  % dont understand this
else
	samplerate = CSC(1).SampleFrequency(1);

end

% clean up memory
clear tempSampleFrequency
clear tempSamples
clear tempTimestamp
