function EF_SpikeSort(fl, DIG, flt, sub, per)

if nargin<5
    per=1;
end
if nargin<4
    sub=[];
elseif ~isempty(sub)
    load (sub)
    sub=raw;
end
if nargin <3
    flt=1;
end
if nargin<2
    DIG=[];
end
s=length(fl);
switch s
    case 1
        sp=[1 1];
        psp={1};
    case 2
        sp=[2 1];
        psp={1; 2};
    case 3
        sp=[2 2];
        psp={[1 2]; 3; 4};
    case 4
        sp=[2 3];
        psp={[1 2 3]; 4; 5; 6};
    otherwise
        warndlg('Too many files!', 'Error','modal');
        return;
end
for i=1:s
    if fl{i}(end-3)=='.'
        fl{i}=fl(1:end-4);
    end
    
    load([fl{i} '.mat']);

    if ~exist('raw','var') || ~exist ('info','var')
        warndlg ('Required data is not present. Please check the file name','ERROR');
        return;
    end
    if sub
        dat=double(raw-sub);
    else
        dat=double(raw);
    end
    switch (flt)
        case 1
            datx=MC_spHighPass(dat,info.Fs);
            dat=datx;
        case 2
            datx=dat-MC_LowPass (dat,info.Fs,3000);
            dat=MC_spHighPass(dat,info.Fs);
        otherwise
        
    end
    
    eval(['dat' num2str(i) '=dat;']);
    eval(['datx' num2str(i) '=dat;']);
end

ff=figure;
subplot(2,1,1);
dd=datx1(round(info.Fs*10):round(info.Fs*12));
plot (dd);
title ('Raw, filtered');
subplot (2,1,2);
plot (dd.*dd);
title ('RMS');
a=questdlg ('Choose method for spike detection','Choose method','Raw','RMS','Raw');
if strcmp (a,'RMS')
    for i=1:s
        dd=eval(['datx' num2str(i)]);
        dd=dd.*dd;
        eval(['datx' num2str(i) '=-dd;']);
    end
end
close (ff);
clear dat datx raw dd;
fls=[fl{1} '_spk.mat'];
ch=ones(1,s);
scrsz = get(0,'ScreenSize');
f=figure('Position',[round(scrsz(3)/10) round(scrsz(4)/10) round(8*scrsz(3)/10) round(8*scrsz(4)/10)]);

rep=1;
l=length(datx1(round(info.Fs*5):round(info.Fs*10)));
yn=min(datx1(round(info.Fs*5):round(info.Fs*10)));
yx=max(datx1(round(info.Fs*5):round(info.Fs*10)));
for i=1:s
    subplot(sp(1),sp(2),psp{i});
    eval(['dat=datx' num2str(i) ';']);
    plot (dat(round(info.Fs*5):round(info.Fs*10)),'k');
    axis ([1 l yn*2 yx*2]);
end
pps=per;
while rep
    yn=min(datx1(round(info.Fs*pps+info.Fs*5):round(info.Fs*pps+info.Fs*10)));
    yx=max(datx1(round(info.Fs*pps+info.Fs*5):round(info.Fs*pps+info.Fs*10)));
    for i=1:s
        subplot(sp(1),sp(2),psp{i});
        eval(['dat=datx' num2str(i) ';']);
        plot (dat(round(info.Fs*pps+info.Fs*5):round(info.Fs*pps+info.Fs*10)),'k');
        axis ([1 l yn*2 yx*2]);
    end
    subplot(sp(1),sp(2),psp{1});
    hold off
    plot (datx1(round(info.Fs*pps+info.Fs*5):round(info.Fs*pps+info.Fs*10)));
    axis ([1 l yn*2 yx*2]);
    hold on;
    [x oy1]=ginput(1);
    line([1 info.Fs*5], [oy1 oy1],'Color','r');
    [x oy2]=ginput(1);
    line([1 info.Fs*5], [oy2 oy2],'Color','r');
    [x ts]=ginput(1);
    line([1 info.Fs*5], [ts ts],'Color','g');
    a=questdlg ('Accept values?', 'Question...','Yes', 'No', 'Next Epoch', 'No');
    if strcmp(a,'Yes')
        rep=0;
    elseif strcmp(a,'Next Epoch')
        pps=pps+60;
    end
end
close (f);
outl=[min(oy1,oy2) max(oy1,oy2)];

[wav spkTimes] =cutSpikes (datx1,info.Fs,ts,outl, 0:(1000/info.Fs):(length(dat1)*1000/info.Fs));
% spkTimes=round(1000*spkTimes)/1000;
a=round(spkTimes/1000*info.Fs);
for i=1:length(spkTimes)
    wav(:,i)=dat1(a(i)-25:a(i)+47);
end
if s>1
%     [x a ~]=intersect(round(1000*(0:(1000/info.Fs):(length(dat1)*1000/info.Fs)))/1000,spkTimes);
    for i=1:length(spkTimes)
        for j=2:s
            dat=eval(['dat' num2str(j)]);
             wav(74+73*(j-2):73+73*(j-1),i)=dat(a(i)-25:a(i)+47);
        end
    end
end
clear dat*
[SpikeTimes,SpikeShapes,Grades]=EF_spSort(wav, spkTimes, DIG);
save (fls,'SpikeTimes','SpikeShapes','Grades');
return;
end

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

end