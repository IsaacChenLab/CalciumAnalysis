function [more,startSeg,varargout]=MC_funcElectrodesSeg(fileName,ch,seg,func,narg,varargin)


[dataStream,totalMs,Hz]=MC_openFile(fileName);

s1=['(rawData,xtime,Hz'];
for i=1:length(varargin)
    s1=[s1 ',' varargin{i}];
end
s1=[func s ')'];
nout = nargout-2;
s2='[';
for i=1:nout,
    s2=[s2 'varargout{' num2str(i) '} '];
end
s2=[s2 ']'];
s=[s1 s2];

periods=floor(totalMS/seg);
last=mod(totalMS,seg);
for i=1:periods
    % get raw data
    [rawData,xtime]=MC_getElectrodes(fileName,ch, Start+[0 seg]+(i-1)*seg );
    % cut spikes by threshold
    eval(s);
end

return;