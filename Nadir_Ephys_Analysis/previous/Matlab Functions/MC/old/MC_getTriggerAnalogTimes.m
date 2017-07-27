function [anaTimes anaInxs analogData xtime]=MC_getTriggerAnalogTimes(fileName,analogCh,Time,threshold,posNegRmsDc)
% [anaTimes anaInxs analogData]=MC_getTriggerAnalogTimes(fileName,analogCh,Time,threshold,)
% triggers an analog channel and returns the time of triggers in ms, the
% indexes in original sampling and the analog data itself.
% threshold is standard deviations if its value is between 0 and 10,
% otherwise it is the actual value.
% posNegRmsDc is either 'p' for positive crossing, 'n' for negative crossing,
% both can have an 'r' for rms calculation, and can have 'd' for
% subtracting the mean before everything
% default is no rms, positive crossing of 3 standard deviations, and dc
% subtracted

if nargin<3
    Time=nan;
end
if nargin<2,
    analogCh=1;
end
if nargin<4,
    threshold=5;
end
if nargin<5
    rms='no';
    slope='p';
else
    slope='p';
    rms='no';
    dc='yes';
    if ~isnan(findstr(posNegRmsDc,'r'))
        rms='yes';
    end
    if ~isnan(findstr(posNegRmsDc,'n'))
        slope='n';
    end
    if ~isnan(findstr(posNegRmsDc,'d'))
        dc='no';
    end    
end

for i=1:length(analogCh)
    
    [anaData,xtime,Hz]=MC_getAnalog(fileName,analogCh,Time);
    analogData(:,i)=anaData;
    
    if strcmp(dc,'yes')
        anaData=anaData-mean(cData);
    end

    if strcmp(rms,'yes')
        % root mean square on 0.2 ms
        ns= Hz /1000 / 5;
        rms1=anaData.^2;
        rms1=rms1(1:ns*floor(length(rms1)/ns));
        rms2=reshape(rms1,ns,length(rms1)/ns);
        rms3=sqrt(mean(rms2)');
        clear rms1 rms2;
        cData=rms3;
        clear rms3;
    else
        cData=anaData;
        clear anaData;
    end
    
    if threshold > 0 & threshold <=15,
        threshold = mean(cData)+threshold*std(cData);
    end
    
    x=runs(cData >= threshold);
    if strcmp(slope,'p'),
        anaInxs{i}=x(1,:);
    else
        anaInxs{i}=x(1,:)+x(2,:)-1;
    end
    
    if strcmp(rms,'yes')
        anaInxs{i}=(anaInxs{i}-1)*5+3;
    end
    anaTimes{i}=xtime(anaInxs{i});
    anaTimes=1e3*Times(inx)
    clear cData;
end

return;

    

    

