function [] = tdt2mat(sigtype,block,largeFileFlag,channels,savewhere,tank)

% [] = tdt2mat(sigtype,block,longFileFlag,channels,savewhere,tank)
% converts data in the TDT data tank into .mat files
% sorted by channels and signal type (lfps, spike times/waveforms etc)
%
% INPUT:
%   sigtype: string to select which signals to convert
%       'lfp','spike' (pre-thresholded), 'raw' (raw spiking trace),
%       'stim' (photodiode), 'speed' (running speed),'all'(default)
%   block: string with the name of the recording (default:
%       name of folder in current directory).
%   longFileFlag = flag to divide long files in subparts, def. = 0
%   channels: array containing channel numbers, default: 1:32
%   savewhere: 'remote' for pc over the network, 'local' for this one,
%       default = remote
%   tank: path to data dank, default C:\Data
%
% Lucas Pinto, July 2010, lucaspinto@berkeley.edu, partly adapted from
% Aaron Kaluszka's code


if nargin == 0
    sigtype = 'all';
    block = cd;
    block = block(10:end);
    largeFileFlag = 0;
    channels = 1:32;
    savewhere = 'local';
    tank = 'I:\EPhys\EPhys';
elseif nargin == 1
    block = cd;
    block = block(9:end);
    largeFileFlag = 0;
    channels = 1:32;
    savewhere = 'local';
    tank = 'I:\EPhys\EPhys';
elseif nargin == 2
    largeFileFlag = 0;
    channels = 1:32;
    savewhere = 'local';
    tank = 'I:\EPhys\EPhys';
elseif nargin == 3
    channels = 1:32;
    savewhere = 'local';
    tank = 'I:\EPhys\EPhys';
elseif nargin == 4
    savewhere = 'local';
    tank = 'I:\EPhys\EPhys';
elseif nargin == 5
    tank = 'I:\EPhys\EPhys';
end

% save mat files in my desktop and create folders if necessary
if strcmpi(savewhere,'remote')
    currdir = '\\lucasdesk\Data';
elseif strcmpi(savewhere,'local')
    f=find (tank=='\');
    currdir = [tank(1:f(end)) 'To Matlab'];
end
cd(currdir);
% d = date;
% [y,m,d] = datevec(d);
% y = num2str(y); y = y(3:4);
% recday = sprintf('%s%02d%02d',y,m,d);
% filename = dir(recday);
% if size(filename,1) == 0
%     mkdir(recday);
% end
% currdir = pwd;
% cd([currdir '\' recday]);
filename = dir(block);
if size(filename,1) == 0
    mkdir(block);
end
currdir = pwd;
cd([currdir '\' block]);

f = figure('visible','off');
info.recName = block;
TT = actxcontrol('TTank.X');
TT.ConnectServer('Local', 'MyClient');
TT.OpenTank(tank, 'R');
TT.SelectBlock(block);
TT.SetGlobals('WavesMemLimit = 2147483647');

% 1) get timestamps for laser and visual stim
% filename = dir('laser_times.mat'); 
% if size(filename,1) == 0
%     tag = 'lasr/';
%     channel = 0;
%     [laser_times] = tdtFromTank(tag, channel, info,TT);
% %     if isnan(laser_times)
% %         try
% %         [startTime,endTime,laser_times] = ascii2mat_times(block);
% %         catch
% %             disp('laser_times is NaN, please export it manually')
% %         end
% %     end
%     save laser_times laser_times
% end

% filename = dir('startTime.mat');
% if size(filename,1) == 0
%     tag = 'stat/';
%     channel = 0;
%     [startTime] = tdtFromTank(tag, channel, info,TT);
% %     if isnan(startTime)
% %         [startTime] = ascii2mat_times(block); % insert code for exporting here
% %     end
%     save startTime startTime
% else load startTime
% end
% 
% filename = dir('endTime.mat');
% if size(filename,1) == 0
%     tag = 'endt/';
%     channel = 0;
%     [endTime] = tdtFromTank(tag, channel, info,TT);
% %     if isnan(startTime)
% %         [startTime,endTime] = ascii2mat_times(block); % insert code for exporting here
% %     end
%     save endTime endTime
% else load endTime
% end

% if Large file, set times for dividing signal vectors
if largeFileFlag
    if isnan(endTime)
        if strcmpi(block(8:10),'nat')
            a = 19;
        elseif strcmpi(block(8:10),'ori')
            a = 23;
        elseif strcmpi(block(8:10),'con')
            a = 28;
        else 
            a = 20;
        end
        approxDur = a*60;
    else approxDur = endTime/24414.06 + 2; % determine approximate recording duration in sec
    end
    nchunks = 20;
    chunk = round(approxDur/nchunks); % divide in nchunks time chunks
    onesamp = 1/24414.06;
    filetimes = zeros(nchunks,2); % start and end times (cols) of each chunk (rows)
    for k = 1:nchunks
        if k == 1
            filetimes(k,:) = [0 chunk];
        elseif k == nchunks
            filetimes(k,:) = [(k-1)*chunk+onesamp 0];
        else filetimes(k,:) = [(k-1)*chunk+onesamp k*chunk];
        end
    end
    info.nchunks = nchunks;
else filetimes = [0 0];
    info.nchunks = 1;
end

% 2) convert and save lfps
if strcmpi(sigtype,'all') || strcmpi(sigtype,'lfps')
    fprintf('converting LFP files')
    tag = 'lfps';
    for ii = channels
        fprintf('.')
        info.channelid = ii;
        if largeFileFlag
            for jj = 1:4
                [lfp, info] = tdtFromTank(tag,ii,info,TT,filetimes(jj,:));
                save(sprintf('%s%02d_part%01d',tag,ii,jj),'lfp','info');
                clear lfp
            end
        else
            [lfp, info] = tdtFromTank(tag,ii,info,TT,filetimes);
            save(sprintf('%s%02d',tag,ii),'lfp','info');
            clear lfp
        end
    end
    fprintf('\n')
end

if strcmpi(sigtype,'digi')
    tag='digi';
    for ii = channels
        fprintf('.')
        info.channelid = ii;
        [digi, info] = tdtFromTank(tag,ii,info,TT,filetimes);
        save(sprintf('%s%02d',tag,ii),'digi','info');
        clear digi;
    end
end

% 3) convert and save raw spiking data
if strcmpi(sigtype,'all') || strcmpi(sigtype,'raws')
    fprintf('converting raw spiking files')
    tag = 'raws';
    for ii = channels
        fprintf('.')
        info.channelid = ii;
        if largeFileFlag
            for jj = 1:size(filetimes,1)
                [raw, info] = tdtFromTank(tag,ii,info,TT,filetimes(jj,:));
                eval(['raw_part' num2str(jj) '=raw;'])
                clear raw
                if jj == 1
                    save(sprintf('%s%02d',tag,ii),sprintf('raw_part%01d',jj),'info');
                else
                    save(sprintf('%s%02d',tag,ii),sprintf('raw_part%01d',jj),'info','-append');
                end
                eval(['clear raw_part' num2str(jj)])
            end
        else
            [raw, info] = tdtFromTank(tag,ii,info,TT,filetimes);
            save(sprintf('%s%02d',tag,ii),'raw','info');
            clear raw
        end
    end
    fprintf('\n')
end

if strcmpi(sigtype,'all') || strcmpi(sigtype,'rawn')
    fprintf('converting raw spiking files')
    tag = 'rawn';
    for ii = channels
        fprintf('.')
        info.channelid = ii;
        if largeFileFlag
            for jj = 1:size(filetimes,1)
                [raw, info] = tdtFromTank(tag,ii,info,TT,filetimes(jj,:));
                eval(['raw_part' num2str(jj) '=raw;'])
                clear raw
                if jj == 1
                    save(sprintf('%s%02d',tag,ii),sprintf('raw_part%01d',jj),'info');
                else
                    save(sprintf('%s%02d',tag,ii),sprintf('raw_part%01d',jj),'info','-append');
                end
                eval(['clear raw_part' num2str(jj)])
            end
        else
            [raw, info] = tdtFromTank(tag,ii,info,TT,filetimes);
            save(sprintf('%s%02d',tag,ii),'raw','info');
            clear raw
        end
    end
    fprintf('\n')
end

% % 4) convert and save spike time and waveform data
% if strcmpi(sigtype,'all') || strcmpi(sigtype,'spike')
%     fprintf('\nconverting spike time and waveform files')
%    tag = 'spks';
%     for ii = channels
%         fprintf('.')
%         info.channelid = ii;
%         [data, info] = tdtFromTank(tag,ii,info,TT,filetimes);
%         spkwav = data.spkwav;
%         spktimes = data.spktimes;
%         save(sprintf('%s%02d',tag,ii),'spkwav','spktimes','info');
%         clear data spkwav spktimes
%     end
% fprintf('\n')
% end

% 5) convert and save photodiode data
if strcmpi(sigtype,'all') || strcmpi(sigtype,'stim') || strcmpi(sigtype,'ctrl')
    fprintf('\nconverting photodiode data...\n')
    tag = 'ctrl';
    channel = 1;
    info.channelid = channel;
    if largeFileFlag
        for jj = 1:size(filetimes,1)
            [stim, info] = tdtFromTank(tag,channel,info,TT,filetimes(jj,:));
            eval(['stim_part' num2str(jj) '=stim;'])
            clear stim
            if jj == 1
                save('stim',sprintf('stim_part%01d',jj),'info');
            else
                save('stim',sprintf('stim_part%01d',jj),'info','-append');
            end
            eval(['clear stim_part' num2str(jj)])
        end
    else
        [stim, info] = tdtFromTank(tag, channel, info,TT,filetimes);
        save('stim','stim','info');
        clear stim
    end
end

% 6) convert and save running speed data
if strcmpi(sigtype,'all') || strcmpi(sigtype,'ctrl')
    fprintf('\nconverting speed data...\n')
    tag = 'ctrl';
    channel = 2;
    info.channelid = channel;
    [speed, info] = tdtFromTank(tag, channel, info,TT,filetimes);
    save('speed','speed','info');
    clear speed
end

% for siyu's thing
if strcmpi(sigtype,'PDec')
    fprintf('converting PDec')
    tag = sigtype;
    for ii = channels
        fprintf('.')
        info.channelid = ii;
        [PDec, info] = tdtFromTank(tag,ii,info,TT,filetimes);
        save(sprintf('%s%02d',tag,ii),'PDec','info');
        clear PDec
    end
    fprintf('\n')
end

if strcmpi(sigtype,'BeMt')
    fprintf('converting BeMt')
    tag = sigtype;
    for ii = channels
        fprintf('.')
        info.channelid = ii;
        [BeMt, info] = tdtFromTank(tag,ii,info,TT);
        save(sprintf('%s%02d',tag,ii),'BeMt','info');
        clear BeMt
    end
    fprintf('\n')
end

TT.CloseTank;
TT.ReleaseServer;
close(f);
end

%% nested function for retrieving data and info from tank
function [data, info] = tdtFromTank(tag, channel, info, TT, filetimes)

if nargin < 5
    filetimes = [0 0];
end

TT.SetGlobalV('Channel',channel);
TT.SetGlobalV('T1',filetimes(1));
TT.SetGlobalV('T2',filetimes(2));
N = TT.ReadEventsV(9999999, tag, channel, 0, filetimes(1),filetimes(2), 'ALL'); %retrieve events
waves = [TT.ParseEvV(0, N)]'; % get data matrix (snipets)
info.ntpoints = size(waves,2); % n points per snipet
info.Fs  = TT.ParseEvInfoV(1, 1, 9); % sampling frequency

if strcmpi(tag,'spks')
    data.spktimes = TT.ParseEvInfoV(0, N, 6);
    data.spkwav = waves;
elseif strcmpi(tag,'lasr/') || strcmpi(tag,'stat/') || strcmpi(tag,'endt/')
    data = TT.ParseEvInfoV(0, N, 7);
else
    data = TT.ReadWavesV(tag); % get data as single vector
end

end