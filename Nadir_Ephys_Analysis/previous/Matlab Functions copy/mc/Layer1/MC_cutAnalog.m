function [anaTimes]=MC_cutAnalog(rawData,Times,Fs,threshold,prd)
% function [anaTimes]=MC_cutAnalog(rawData,Times,threshold,prd)
% searches for threshold crossings in an analog channel and returns the
% time of triggers in ms.
% Times is the original times in seconds corresponding to the data in
% threshold is standard deviations if its value is between 0 and 15,
% otherwise it is the actual value.
% prd is either 'p' for positive crossing, 'n' for negative crossing,
% both can have an 'r' for rms calculation.
% default is no rms, positive crossing of 5 standard deviations

if nargin<4,
    threshold=5;
end
if nargin<5
    rms='no';
    slope='p';
else
    slope='p';
    rms='no';
    if ~isnan(findstr(prd,'r'))
        rms='yes';
    end
    if ~isnan(findstr(prd,'n'))
        slope='n';
    end
end


if strcmp(rms,'yes')
    [cData,rmsTime,rmsSamples]=MC_rms(rawData,Times,Fs);
else
    cData=rawData;
end
clear rawData;

if threshold > 0 & threshold <=15,
    threshold = mean(cData)+threshold*std(cData);
end

x=runs(cData >= threshold);
if strcmp(slope,'p'),
    Inx=x(1,:);
else
    Inx=x(1,:)+x(2,:)-1;
end

if strcmp(rms,'yes')
    Inx=(Inx-1)*rmsSamples+round(rmsSamples/2); % middle of the 0.2ms (5 samples)
end
% get rid of the second volley if its within 10 ms
SECOND_VOLLEY=10; % 
refractory_samples=SECOND_VOLLEY*Fs/1e3;
df=diff(Inx); x=find(df<refractory_samples); Inx(x+1)=[];

anaTimes=Times(Inx);
clear cData;

return;

    

    

