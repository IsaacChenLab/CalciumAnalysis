function [Inxs]=MC_spCluster(waveForms)
% function [Clusters_SpikeTimes,Clusters_SpikeShapes,Clusters_Mat,Clusters_OriginalIdx]=MC_spikeSort(waveForms,spikeTimes,p1,p2)
% performs spike sorting using unsupervised clustering. first cleans the
% data from outliers, realigns the wave forms according to their minimum,
% then performs unsupervised clustering, and allows merging ot deleting
% groups.
% inputs:
% waveForms is a matrix of -2 to +3ms shapes of spikes sampled at 25e3 Hz
% (samples * #_of_spikes)
% spikeTimes is a vector of all spikeTimes (threshold crossings as detected
% by either MC_cutSpikes or MC_Rack)
% p1 and p2 can be either a number for a pre-determined number of clusters,
% or/and a method for pre-processing for clustering.
% if number of clusters is not supplied, then the most appropriate one,
% based on Gap analysis (Tibshirani and Hastie), is used.
% methodss can be either: 'pca' - takes the ceofficients that extarct 90%
% of the variance; 'shapes' (Default) - origianl wave forms; 'similarity' - creates a
% matrix of similarities between all shapes using the cosine of the angle betwen the
% vectors (as in Sejnowski); 'corrmat' - similarities based on pearson
% correlation; 'spectral' - spectral clustering based on Ng. Weiss. Jordan.
% (NIPS,2001)
% outputs:
% Clusters_SpikeTimes: cell array of spike times for each group
% Clusters_SpikeShapes: cell array of spike wave forms for each group
% Clsuters_Mat: a sparse matrix where each columns represent the occurences
% of a group, 1/0.
% Clusters_OriginalIdx: cell array with indices to the input

if isempty(waveForms)
    error('no spikes provided');
end

Distance='sqEuclidean';
method='decide';

[nSamples,nSpikes]=size(waveForms);
mat=waveForms';
if strcmp(method,'decide')
    if nSpikes<nSamples,
        method='similarity';
    else
        method='pca';
    end
end

method=questdlg('method?','','pca','spectral','similarity','pca');

nClust=MC_getNumClust();


% prepare matrix for clustering
if strcmp(method,'pca'),
    s=inputdlg({'# of coefficients?'},'',1,{'5'});
    NumberCoef=str2num(s{1});
    [mat_c,pcs,expect,lat,coef]=MC_spPCA(mat,NumberCoef);
elseif strcmp(method,'corrmat');
    mat_t=corr(mat','rows','complete');
    R=2*sum(sum(tril(mat_t,-1))) / (size(mat_t,1)*(size(mat_t,1)-1)) ;
    %input reshape
    mmc=mean(mean(mat_t));
    mat_c= 1./ ( 1 + exp(-1*( (mat_t-mmc) / 0.17)) );
elseif strcmp(method,'dot')
    mat_c=mat*mat';
elseif strcmp(method,'similarity')
    %similarity matrix
    mat_t=MC_similarityMatrix(mat);
    %R=2*sum(sum(tril(mat_t,-1))) / (size(mat_t,1)*(size(mat_t,1)-1)) ;
    %input reshape
    mmc=mean(mean(mat_t));
    mat_c= 1./ ( 1 + exp(-1*( (mat_t-mmc) / 0.17)) );
elseif strcmp(method,'shapes')
    mat_c=mat;
elseif strcmp(method,'spectral')
    sigma=10;
    s=inputdlg({'sigma?'},'',1,{'10'});
    sigma=str2num(s{1});
    mat_c = MC_spSpectralEmbed(nClust,mat,sigma);
else
    error('method not recognized');
end

% clustering
if nClust>1
    [idx_c,cent_c,sumdist_c] = kmeans(mat_c,nClust,'display','notify','replicates',20,'distance',Distance);
else
    idx_c=ones(1,size(mat_c,1));
end

for i=1:nClust,
    Inxs{i}=find(idx_c==i);    
end

return;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nClust]=MC_getNumClust()


nClust=questdlg('Number of clusters?','','2','3','4','2');
nClust=str2num(nClust);
if nClust<2 | nClust>4
    error;
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S=MC_similarityMatrix(A)
    
    na=sqrt(sum(A'.^2));
    NA=na'*na;
    S = A*A';
    S=S ./ NA;

return;
