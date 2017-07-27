function [times,nMat,oMat]=MC_relevantDigTimes(type,file,total_length,K)

if nargin<4
    K=50;
end
if strcmp(type,'all')
    times=nan;
    nMat=nan; oMat=nan;
    return;
end

nMat=[]; oMat=[];

if strncmp(type,'rew',3) | strncmp(type,'cs',2),
  if strncmp(type,'rew',3)
      BIT=1;
  elseif strncmp(type,'cs+',3)
      BIT=2;
  elseif strncmp(type,'cs-',3)
      BIT=3;
  elseif strncmp(type,'css',3)
      BIT=-1;
  else
      error('');
  end
  from=str2num(type(4:9));
  to=str2num(type(10:15));
  [Mat,Times]=MC_loadEvent(file);  
  if ~isempty(Mat)
      if BIT==-1
        x=find(Mat(:,2) | Mat(:,3));
      else
        x=find(Mat(:,BIT));  
      end
      nMat=sparse(size(Mat,1),1);
      for i=1:length(x),
          if x(i)+from <1
              ;
          elseif x(i)+to > length(nMat)
              nMat(x(i)+from : end)=1;
          else
              nMat(x(i)+from : x(i)+to)=1;
          end
      end
  end      
 else
    switch type(1:3)
        case 'all'
          if length(type)>3
              [nMat,oMat]=MC_randSelect(ones(total_length,1),str2num(type(4:9)),K);
          end
        case 'rem'
          [nMat,Times]=MC_loadREM(file);  
          if length(type)>3
              [nMat,oMat]=MC_randSelect(nMat,str2num(type(4:9)),K);
          end
        case 'sws'
          [nMat,Times]=MC_loadSWS(file);  
          if length(type)>3
              [nMat,oMat]=MC_randSelect(nMat,str2num(type(4:9)),K);
          end
    end
end

times=find(nMat);
if isnan(times)
    times=[];
end
if isnan(nMat)
    nMat=[];
end
if isempty(oMat) | isnan(oMat)
    oMat=nMat;
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%
function [fMat,oMat]=MC_randSelect(nMat,segl,K)

if isempty(nMat)
    oMat=[]; fMat=[];
    return;
end

%K=150; % max 100 segments

x=find(nMat);
n=length(x);
t=floor(n/(2*segl));
K=min(t,K);
al=round(linspace(random('unid',segl,1,1),n-2*segl,K));
al=al+random('unid',segl,1,K);
oMat=sparse(size(nMat,1),1);
for i=al
    oMat(x(i):x(i)+segl-1,1)=1;
end

fMat=sparse(oMat&nMat);

return;



