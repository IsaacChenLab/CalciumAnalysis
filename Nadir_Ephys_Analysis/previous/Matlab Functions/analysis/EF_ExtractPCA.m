function [pc, exp]=EF_ExtractPCA (mat)

[i,~]=find(~isnan(mat));
mato=mat(unique(i),:);
matcov=cov(mato);
[pc,~,exp]=pcacov(matcov);