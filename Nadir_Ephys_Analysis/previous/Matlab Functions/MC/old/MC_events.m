function [EventTimes,Mat,Name]=MC_events(file,eventType,p1,Time)
% function [EventTimes,Mat,Name]=MC_events(fileName,eventType,p1,Time)
% extract digital events (up bit) or analog event (threshold crossings) by
% files.
% eventType can be 'dig' or 'ana', and p1 is either the digital bit number
% (can be an array of some), or the analog channel number (can be an array)
% returns a cell array of event times in seconds, a sparse matrix of length
% of file in ms times number of events requested with one where an event
% has occured; and the names of the events.

if nargin<5,
    Time=nan;
end

file=MC_fileNumber(file);

Total_ms=0; EventTimes=[];
% unite events over files
for i=1:length(file)
    fileName=MC_fileName(file(i));
    [dataStream,totalMS,Hz]=MC_openFile(fileName);
    if strncmp(eventType,'dig',3)
        bits=p1;
        nEvent=length(bits);
        [st]=MC_getDigitalTimes(fileName,bits,Time);
        nn='digital';
    elseif strncmp(eventType,'ana',3)
        anaCh=p1;
        nEvent=length(anaCh);
        [rawData,xtime]=MC_getAnalog(fileName,anaCh,Time);
        [st]=MC_cutAnalog(rawData,xtime,Hz,10,'pr');
        clear rawData xtime;
        nn='analog';
    end
    EventTimes=[EventTimes; Total_ms+st];
    Total_ms=Total_ms+totalMS;
end

% create final matrix
Mat=sparse(Total_ms,nEvent);
for i=1:nEvent
    Mat(round(EventTimes),i)=1;
    Name{i}=sprintf('%s_%02d',nn,p1(i));
end
   
return;

    

    
