function [partition, obj] = graclus(G, k, cutType, l, spectral)
%% Graclus Matlab interface
%
% [partition, obj] = graclus(G, k, cutType, l, spectral)
%
% G: Graph adjacency matrix
% k: number of clusters
% cutType: 0 for NCut, 1 for RAssoc (default is NCut)
% l: number of local search steps (default is 0)
% spectral: 1 for spectral clustering at coarsest level, 0
% otherwise (default is not to use spectral clustering)
if nargin < 5
  spectral = 0;
end
if nargin < 4
    l = 0;
end
if nargin < 3
    cutType = 0;
end
if nargin < 2
    disp('Please specify the adjacency matrix and the number of clusters');
    partition = 0;
    obj = 0;
else
    if ~issparse(G)
      G = sparse(G);
    end
    %Graclus only accepts integer-valued edge weights, so round here
    %Note: may need to round more if more digits of accuracy are needed
    G = round(1000*G);
    [partition, obj] = graclus_mex(G,nnz(G),k,cutType,l,spectral);
end