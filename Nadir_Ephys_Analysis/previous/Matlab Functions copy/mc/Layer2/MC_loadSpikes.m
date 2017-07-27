function [SpikeMat,SpikeGrade,SpikeInxs,SpikeTimes,SpikeShapes]=MC_loadSpikes(file,ch,sp)

% sp==-1 chooses the best spike (first by grade then by number of spikes)

global CAT;
global DATE;
global DATA_DIR;

SpikeMat=[];
SpikeGrade=[];
SpikeTimes=[];
SpikeShapes=[];
SpikeInxs=[];

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  
nfl=MC_multiFilesName(file);
fileName=MC_digFileName(fileName);
spikeShapes_file=MC_spShapeFileName(file);

s1=(sprintf('SpikeMat_%d',ch));
s2=(sprintf('SpikeGrades_%d',ch));
s3=(sprintf('SpikeTimes_%d',ch));
s4=(sprintf('SpikeShapes_%d',ch));

s=sprintf('load %s %s %s %s;',fileName,s1,s2,s3);
warning off all;
eval(s);
if nargout==5
    s=sprintf('load %s %s;',spikeShapes_file,s4);
    eval(s);
end
warning on;

if ~exist(s1,'var') | isempty(eval(s1)) | isnan(eval(s1));
    return;
end

SpikeMat=eval(s1);
SpikeGrade=eval(s2);
SpikeTimes=eval(s3);
SpikeInxs=1:size(SpikeMat,2);

if nargout==5
    SpikeShapes=eval(s4);
end

if nargin==3 & sp>0
    x=find(SpikeGrade<=sp);
    SpikeMat=SpikeMat(:,x);
    SpikeGrade=SpikeGrade(x);
    SpikeTimes=SpikeTimes(x);
    SpikeInxs=SpikeInxs(x);
    if nargout==5
        SpikeShapes=SpikeShapes(x);
    end
elseif nargin==3 & length(SpikeGrade)>1 & sp==-1
    [j,i]=sort(SpikeGrade,'ascend');
    x=find(j==j(1));
    if length(x)>1
        for k=1:length(x),
            l(k)=length(SpikeTimes{i(k)});
        end
        [jj,ii]=sort(l,'descend');
        sp=i(ii(1));
    else
        sp=i(1);
    end
    SpikeMat=SpikeMat(:,sp);
    SpikeGrade=SpikeGrade(sp);
    SpikeTimes=SpikeTimes{sp};
    SpikeInxs=SpikeInxs(sp);
    if nargout==5
        SpikeShapes=SpikeShapes(sp);
    end

end
    
return;