function [waveForms,spikeTimes,rms,threshold,outliers]=...
    MC_cutSpikes(rawData,Times,digData,Fs,rms,threshold,outliers,whatevents,evinterv)
% function [waveForms,spikeTimes]=MC_cutSpikes(rawData,Times,Fs,Threshold)
% filters the rawData with a high pass filter, finds thresholds crossings,
% and extracts -1 till +2 ms as the spikes waveforms.
% inputs:
% rawData is a vector of voltages
% Times is the corresponding time in second
% Fs assumed to be 25e3, unless supplied
% procedure:
% 1. filters high pass
% 2. shows data and asks wether RMS or not should be used
% 3. asks for outliers threshold i.e. a y-axis value that if corssed, the
% 'fake spike' will not counted.
% 4. asks for threshold, i.e. a y-axis value that if corossed the event
% will be counted as spike.
% 5. extract the spikes wave forms, -2ms to +3ms.
% 6. realign them on minimum voltage and takes -1ms to +2ms

waveForms=[]; spikeTimes=[];
if nargin<4,
    Fs=25e3;
end
SECOND_VOLLEY=2; % refractory of 1 ms
refractory_samples=SECOND_VOLLEY*Fs/1e3;
MINUS=2; PLUS=3; % minus plus times taken from each spike
minus_samples=MINUS*Fs/1e3;
plus_samples=PLUS*Fs/1e3;

if nargin<5
% choose thresholds and rms 
    [rms,threshold,outliers,whatevents,evinterv]=MC_CutSpikesGetParam(rawData,Times,digData,Fs);
end

fltData=MC_spHighPass(rawData,Fs);
clear rawData;

%% remove the intervals for events
[fltData]=MC_spRemoveAroundEvents(fltData,digData,Times,whatevents,evinterv);

% find spikes
if rms
    [rmsData,rmsTime,rmsSamples]=MC_spRMS(fltData,Times,Fs);
    [inx]=find(rmsData(:,1)>threshold);
    inx=MC_runs(inx'); % find the beginning of each run
    inx=inx(1,:);
    % get rid of the second volley if its within 2 ms, and edges
    refractory_samples=refractory_samples/rmsSamples;
    df=diff(inx); x=find(df<refractory_samples); inx(x+1)=[];
    inx=inx(inx-minus_samples/rmsSamples >0 & inx+plus_samples/rmsSamples<length(rmsData));    
    if ~isempty(inx)
        % create the matching wave forms for spikes
        v=repmat([-minus_samples/rmsSamples:plus_samples/rmsSamples]',1,length(inx));
        waveInx=v+repmat(inx,length(-minus_samples/rmsSamples:plus_samples/rmsSamples),1);
        waveForms=rmsData(waveInx);
        clear rmsData;
        % outliers
        inx=inx(find((max(waveForms)<outliers)));
        % return to original
        inx=(inx-1)*rmsSamples+round(rmsSamples/2); % middle of the 0.2ms (5 samples)
    end
else
    [inx]=find(fltData(:,1)<threshold);
    inx=MC_runs(inx'); % find the beginning of each run
    inx=inx(1,:);
    % get rid of the second volley if its within 2 ms, and edges
    df=diff(inx); x=find(df<refractory_samples); inx(x+1)=[];
    inx=inx(inx-minus_samples >0 & inx+plus_samples<length(fltData)); 
    if ~isempty(inx)
        % create the matching wave forms for spikes
        v=repmat([-minus_samples:plus_samples]',1,length(inx));
        waveInx=v+repmat(inx,length(-minus_samples:plus_samples),1);
        waveForms=fltData(waveInx);
        % outliers
        inx=inx( find( (min(waveForms)>outliers(1)) & (max(waveForms)<outliers(2)) ) );
    end
end    

if ~isempty(inx)
    v=repmat([-minus_samples:plus_samples]',1,length(inx));
    waveInx=v+repmat(inx,length(-minus_samples:plus_samples),1);
    waveForms=fltData(waveInx);

    spikeTimes=Times(inx); % in ms

    % realign shapes by minimum
    [waveForms,inx,shifts]=MC_spRealign(waveForms,Fs);
    spikeTimes=spikeTimes(inx)+[1e3*shifts/Fs]';
end

return;

