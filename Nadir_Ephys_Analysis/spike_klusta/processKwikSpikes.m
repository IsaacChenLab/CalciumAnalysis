function mydv = processKwikSpikes(filebase,numdvs,group,tvectors, datablock,samplerate)

% returns a struct will lots of good info about clusters at each dv

begin_win = 15;
end_win = 16;

for j=1:numdvs
    disp('DV:');
    disp(j);
    if numdvs==1
        [sp st cl feat gr]=getKwikSpikeData(filebase,tvectors(j).tvector);
    else
        [sp st cl feat gr]=getKwikSpikeData([filebase num2str(j)],tvectors(j).tvector);
    end
    %wv=getAllWaves([filebase num2str(j) '.kwx']); % outdated
    wv=getWaves(sp,datablock,begin_win,end_win);
    unt = getAllUnits(sp,st,cl,feat,wv,gr,group,samplerate);
dv(j)=struct('spikes',sp,'stamps',st,'clusters',cl,'features',feat,'groups',gr,'waves',wv,'units',unt);
end

mydv=dv;

end