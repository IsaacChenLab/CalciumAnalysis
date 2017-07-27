function MC_convertToBin(files,channels,ELEC_CHUNK)

global  DATA_DIR;
global CAT;
cdc;
datapath=[DATA_DIR '\' CAT '\IgorData\'];

if nargin<3
    ELEC_CHUNK=4;
end
MAX_CHANNELS=44;
CHUNKS=5000; % chunks of 5000 ms

filelist=dir([CAT '*.mcd']);
if nargin==0 | files==nan
    files=[1:length(filelist)];
    filenames_supplied=0;
else
    filenames_supplied=1;
end

for ifl=files,
    clear dataStream totalSeconds Sampling_Hz flname elecData;
    if filenames_supplied
        flname=[CAT sprintf('%03d',ifl) '.mcd'];
    else
        flname=filelist(ifl).name;
    end
    mark=(findstr(flname,CAT)+length(CAT));
    q=findstr(flname,'.')-mark;
    if q==0 %
        nfl=0;
    else
        nfl=str2num(flname(mark:mark+2));
    end
    
    [dataStream,totalSeconds,Sampling_Hz]=MC_openFile(flname);
    
    datastrms=getfield(dataStream,'StreamNames');
    datainfos=getfield(dataStream,'StreamInfo');
    nch=getfield(dataStream,'NChannels2');
    
    %%%%%%%%%%%%%%%%
    % electrodes raw data 
    %%%%%%%%%%%%%%%%%%%%%%%%%
    streamName='Electrode Raw Data';
    i=strmatch(streamName,datastrms,'exact');
    cnch=min(nch(i),MAX_CHANNELS);
    if nargin<2
        channels=1:cnch;
    end
    elecData=[];
    totalMS=(totalSeconds)*1000;
    Start=0;
    periods=floor(totalMS/CHUNKS);
    last=mod(totalMS,CHUNKS);
    for i=1:floor(length(channels)/ELEC_CHUNK),
        channelNumbers=channels([1 : ELEC_CHUNK] + (i-1)*ELEC_CHUNK)
        outfilename=sprintf('%s%s%03d_ch%d-%d.bin',datapath,CAT,nfl,channelNumbers(1),channelNumbers(end));
        fid=fopen(outfilename,'w');
        for i=1:periods,
            ce = nextdata(dataStream,'streamname',streamName,'originorder','on','startend',Start+[0 CHUNKS]+(i-1)*CHUNKS);
            [vdata, tvals] = ad2muvolt(dataStream, ce.data(channelNumbers,:), streamName);
            fwrite(fid,vdata','float32');
        end
        ce = nextdata(dataStream,'streamname',streamName,'originorder','on','startend',Start+[periods*CHUNKS periods*CHUNKS+last]);
        [vdata, tvals] = ad2muvolt(dataStream, ce.data(channelNumbers,:), streamName);
        fwrite(fid,vdata','float32');
        fclose(fid); 
    end
    
    % triggers
    
    % ekg
    
    % spikes
    

end

    
return;
