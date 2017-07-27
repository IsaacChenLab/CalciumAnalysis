function [mat_c,pc,exp,lat,coef]=MC_spPCA(mat,NumberCoef)

% global DATA_DIR;

pc=[];lat=[];
if exist('PCA_spike_shapes.mat','file')
    load('PCA_spike_shapes.mat');
end
if nargin<2
    cms=cumsum(exp); 
    x=find(cms>90);
    NumberCoef=max(x(1),5);
end
if size(mat,2)~=size(pc,2)
    pts=size(pc,2);
    n=floor(size(mat,2)/pts);
    mat_c=zeros(size(mat,1),n*NumberCoef);
    for i=1:n
        cf=mat(:,1+(i-1)*pts:i*pts)*pc;
        mat_c(:,1+(i-1)*NumberCoef:i*NumberCoef)=cf(:,1:NumberCoef);
    end
    coef=cf;
else
    coef=mat*pc;
    mat_c=coef(:,1:NumberCoef);
end

return;
