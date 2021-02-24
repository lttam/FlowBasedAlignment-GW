function [K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = figtreeKCenterClustering(d, N, x, K)
%
% (THIRD-PARTY TOOLBOX)
%     Gonzalez's farthest-point clustering algorithm.
% See the supplementary material for more information.
%% Input
%
%    * d --> data dimensionality.
%    * N --> number of source points.
%    * x --> d x N matrix of N source points in d dimensions 
%    * kMax --> maximum number of clusters.
%
%% Ouput
%
%    * K --> actual number of clusters (less than kMax if duplicate pts exist)
%    * rx --> maximum radius of the clusters.
%    * clusterIndex --> vector of length N where the i th element is the 
%                cluster number to which the i th point belongs. 
%                ClusterIndex[i] varies between 0 to K-1.
%    * clusterCenters --> d x K matrix of K cluster centers 
%    * numPoints --> number of points in each cluster.
%    * clusterRadii --> radius of each cluster.
%
