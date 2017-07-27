function [bitTimes,Mat,FileTimes]=MC_evGet(file,Time)

if nargin<2,
    Time=nan;
end

file=MC_fileNumber(file);

bitTimes=cell(1,16);

Total_ms=0; EventTimes=[];
% unite events over files
for i=1:length(file)
    file(i)
    fileName=MC_fileName(file(i));
    [st,xtime,Hz,totalMS]=MC_loadDigital(file(i),Time);
    [b] =  MC_dec2bin(st,xtime);
    for i=1:16,
        bitTimes{i}=[bitTimes{i}; b{i}+Total_ms];
    end
    Total_ms=Total_ms+totalMS;
    FileTimes(i)=Total_ms;
end

% create final matrix
Mat=sparse(Total_ms,16);
for i=1:16
    Mat(round(bitTimes{i}),i)=1;
end
   
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b] =  MC_dec2bin(decvec,time)

%aa=intersect(find(diff(decvec))+1,find(decvec)); % onsets of different decimal numbers
L=length(decvec);

xt=find(decvec);
if isempty(xt)
    b=cell(1,16);
    return;
end
decvec=decvec(xt);
xtime=time(xt);
inxs=1:L;
binstr=dec2bin(decvec,16);
binstr=fliplr(binstr);
binvec=sparse(binstr-48);

for i=1:16,
    xx=find(binvec(:,i));
    x=MC_runs(inxs(xt(xx)));
    x=x(1,:);
    b{i}=time(x);
end

return;



    
