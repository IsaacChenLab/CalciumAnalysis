function [rmsData,xtime,rmsSamples]=MC_spRMS(rawData,xtime,Hz,INTEG)
% function [rmsData,rmsSamples]=MC_rms(rawData,xtime,Hz,INTEG)
% retruns the RMS signal and the number of samples in the original signal
% per one sample in the rms signal. integartes over INTEG ms, default is
% 0.2.

if nargin<4
    INTEG=0.2; % over 0.2ms
end

rmsSamples=round(INTEG*Hz/1e3);

for i=1:size(rawData,2)
    % root mean square on 0.2 ms
    rms1=rawData(:,i).^2;
    rms1=rms1(1:rmsSamples*floor(length(rms1)/rmsSamples));
    rms2=reshape(rms1,rmsSamples,length(rms1)/rmsSamples);
    rmsData(:,i)=sqrt(mean(rms2)');
end
xtime=decimate(xtime,rmsSamples);
xtime=xtime(1:size(rmsData,1));

return;




