function plotUnits(dv, myunits, ADBitVolts, autowidth, recordinglength)

%autowidth = 250; in ms
% samplerate in samples per second
autobins = 100;

numFigs = ceil((size(myunits,1)*2) / 12);
disp(numFigs);
figwidth = 6;
figheight = ceil((size(myunits,1)*2) / figwidth);

j=1;
for k=1:numFigs
    if k*12 > size(myunits,1)*2
        numplots = size(myunits,1)*2 - (k-1)*12;
    else
        numplots = 12;
    end
    disp(numplots);
    figure;
for i=1:2:numplots
    ac=spikeAuto(dv(myunits(j,1)).units(myunits(j,2)).stamps,0);
    ac=ac/1000;
    ad=ac(find(abs(ac)<autowidth));
    fr = length(dv(myunits(j,1)).units(myunits(j,2)).spikes) / recordinglength;
    subplot(3,4,i);
    histogram(ad,autobins);
    title({'Autocorrelation';['DV ',num2str(myunits(j,1)),' UNIT ',num2str(myunits(j,2))]; ['FR: ',num2str(fr)]});
    subplot(3,4,i+1);
    plotMeanWave(dv(myunits(j,1)).units(myunits(j,2)).waves,dv(myunits(j,1)).units(myunits(j,2)).maxChannel,ADBitVolts*1000000);
    title({'Waveform w stdev';['DV ',num2str(myunits(j,1)),' UNIT ',num2str(myunits(j,2))]});
    j=j+1;
end
end
    