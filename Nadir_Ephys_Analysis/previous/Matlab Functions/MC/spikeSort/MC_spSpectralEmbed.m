function mat = MC_spSpectralEmbed(k,M,sigma)
% idx - the partition of the nodes to clusters
% k - number of clusters
% M - n*m data matrix (n data points - each point is a m dimesional row vector)

[n,m] = size(M);
%create the distance matrix
for i = 1:n
    tmp_mat = repmat(M(i,:),n-i,1)-M(i+1:n,:);
    dist_mat(i,i+1:n) = diag(tmp_mat*tmp_mat')';
    clear tmp_mat;
    dist_mat(i+1:n,i) = dist_mat(i,i+1:n)';
end
%create the similarity matrix
A = exp(1).^(-dist_mat/(2*sigma^2));
clear dist_mat;
A = A-eye(n);
%create the laplacian
DD = diag(sum(A));
DD2 = diag(sum(A).^-0.5);
L = DD2*(DD-A)*DD2;
clear A;
%find the k smallest eigen vectors 
[V,D]=eig(L);
[eig_vals,inds_k] = sort(diag(D));
X = V(:,inds_k(1:k));
%cluster the nodes
Y = normr(X);
mat=Y;
