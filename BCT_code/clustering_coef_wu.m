function C=clustering_coef_wu(W)
%CLUSTERING_COEF_WU     Clustering coefficient
%
%   C = clustering_coef_wu(W);
%
%   The weighted clustering coefficient is the average "intensity"
%   (geometric mean) of all triangles associated with each node.
%
%   Input:      W,      weighted undirected connection matrix
%                       (all weights must be between 0 and 1)
%
%   Output:     C,      clustering coefficient vector
%
%   Note:   All weights must be between 0 and 1.
%           This may be achieved using the weight_conversion.m function,
%           W_nrm = weight_conversion(W, 'normalize');
%
%   Reference: Onnela et al. (2005) Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2015

%   Modification history:
%   2007: original
%   2015: expanded documentation


K=sum(W~=0,2);              % K is a column vector of the degree of every node 

cyc3=diag((W.^(1/3))^3);    % cyc3 is an n-length vector containing the sum of the geometric means of the weights of all the complete triangles centered on a given node
                            % the diagonal of the cubed matrix gives you the sum of the products of the weights of the edges along all the possible 3-step paths from a node back to itself
                            %(ie a complete triangle). Taking the 3rd root of the matrix before hand turns that into the sum of the geometric means of the weights of the triangles
                                                      
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)

C=cyc3./(K.*(K-1));         %clustering coefficient
