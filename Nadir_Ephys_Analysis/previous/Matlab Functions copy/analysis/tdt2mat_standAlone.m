function tdt2mat_standAlone(block)

rootdir = pwd;
% block = 'Block-2';
mkdir(block)
cd([rootdir '\' block]);
tank = 'C:\Data\100910';
tag = 'PD__';
channels = 1:2;

f = figure('visible','off');
info.recName = block;
TT = actxcontrol('TTank.X');
TT.ConnectServer('Local', 'MyClient');
TT.OpenTank(tank, 'R');
TT.SelectBlock(block);
TT.SetGlobals('WavesMemLimit = 2147483647');

for ii = channels
    [raw, info] = tdtFromTank(tag,ii,info,TT);
    save(sprintf('%s%02d',tag,ii),'raw','info');
    clear raw
end

TT.CloseTank;
TT.ReleaseServer;
close(f);
cd(rootdir)
end

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