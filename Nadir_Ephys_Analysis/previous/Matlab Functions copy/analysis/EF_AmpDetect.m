function [spikes,index]  = EF_AmpDetect(x,Fs)

dat1=MC_spHighPass(double(x),Fs);

scrsz = get(0,'ScreenSize');
f=figure('Position',[round(scrsz(3)/10) round(scrsz(4)/10) round(8*scrsz(3)/10) round(8*scrsz(4)/10)]);

rep=1;
l=length(dat1(round(Fs*5):round(Fs*10)));
pps=1;
while rep
    yn=min(dat1(round(Fs*pps+Fs*5):round(Fs*pps+Fs*10)));
    yx=max(dat1(round(Fs*pps+Fs*5):round(Fs*pps+Fs*10)));
    hold off
    plot (dat1(round(Fs*pps+Fs*5):round(Fs*pps+Fs*10)));
    axis ([1 l yn*2 yx*2]);
    hold on;
    [x oy1]=ginput(1);
    line([1 Fs*5], [oy1 oy1],'Color','r');
    [x oy2]=ginput(1);
    line([1 Fs*5], [oy2 oy2],'Color','r');
    [x ts]=ginput(1);
    line([1 Fs*5], [ts ts],'Color','g');
    a=questdlg ('Accept values?', 'Question...','Yes', 'No', 'Next Epoch', 'No');
    if strcmp(a,'Yes')
        rep=0;
    elseif strcmp(a,'Next Epoch')
        pps=pps+60;
    end
end
close (f);
outl=[min(oy1,oy2) max(oy1,oy2)];

[wav, spkTimes] =cutSpikes (dat1,Fs,ts,outl, 0:(1000/Fs):(length(dat1)*1000/Fs));

spikes = wav';
index = spkTimes;
function [waveForms,spikeTimes]= cutSpikes(fltData,Fs,threshold,outliers, Times)
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

% find spikes
[inx]=find(fltData<threshold);
inx=MC_runs(inx'); % find the beginning of each run
inx=inx(1,:);
% get rid of the second volley if its within 2 ms, and edges
df=diff(inx); x=find(df<refractory_samples); inx(x+1)=[];
inx=inx(inx-minus_samples >0 & inx+plus_samples<length(fltData)); 
if ~isempty(inx)
    % create the matching wave forms for spikes
    v=repmat([-minus_samples:plus_samples]',1,length(inx));
    waveInx=v+repmat(inx,length(-minus_samples:plus_samples),1);
    waveForms=fltData(round(waveInx));
    % outliers
    inx=inx( find( (min(waveForms)>outliers(1)) & (max(waveForms)<outliers(2)) ) );
end

if ~isempty(inx)
    v=repmat([-minus_samples:plus_samples]',1,length(inx));
    waveInx=v+repmat(inx,length(-minus_samples:plus_samples),1);
    waveForms=fltData(round(waveInx));

    spikeTimes=Times(inx); % in ms

    % realign shapes by minimum
    [waveForms,inx,shifts]=MC_spRealign(waveForms,Fs);
    spikeTimes=spikeTimes(inx)+[1e3*shifts/Fs];
end

return;
