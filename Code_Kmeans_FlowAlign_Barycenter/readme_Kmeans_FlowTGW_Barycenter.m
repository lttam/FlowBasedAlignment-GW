%%%%%%%%%%%%%
% GUIDELINE %
%%%%%%%%%%%%%

% Illustration Kmeans clustering with FlowAlign barycenter
% on 60K MNIST samples rotated randomly
% (FlowAlign: sampling aligned-root tree metric --> So, one can apply
% aligned-root FlowAlign to reduce its complexity)

% -------------------------------------
% Step 1: run "preKmeans_MNISTxx60K_clouds.m"
% + Building tree metrics (same structure for each probability measure on
% different spaces) by clustering-based tree metric sampling
% + Having FlowTreeGW representation

% -------------------------------------
% Step 2: run "compute_Kmeans9_MNIST_xx60K_clouds.m"
% + compute Kmeans with FlowTreeGW barycenters

% -------------------------------------
% Step 3: run "evaluate_Kmeans9_MNIST_xx60K_clouds.m"
% + evaluate the performance of Kmeans clustering by F-beta measure
% -------------------------------------


