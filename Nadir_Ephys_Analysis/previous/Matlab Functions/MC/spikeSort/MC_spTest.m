function MC_spTest

fl=0/1/2; % first or more file
a=find(SpikeMat_1(:,1));
k=50; % spike
totalMs= % total ms in file

fileName=MC_fileName(file);
[rawData,xtime]=MC_getElectrodes(fileName,1,nan);
fltData=MC_spHighPass(rawData,25e3);
subplot(1,2,1); 
plot(SpikeShapes_1{1}(:,k)); 
subplot(1,2,2); 
plot(fltData((a(k)-fl*totalMs)*25-25 : (a(k)-fl*totalMs)*25+50));