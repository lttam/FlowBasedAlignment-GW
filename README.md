# FlowBasedAlignment-GW
Code and data for "Flow-based Alignment Approaches for Probability Measures in Different Spaces" -- AISTATS 2021 (Tam Le, Nhat Ho, Makoto Yamada)


% % Some guidelines for the code of "Flow-based alignment approaches for probability measures in different spaces"
% % 
% % (a) Dataset
% % For an example: twitter_org.mat : TWITTER dataset (documents with word embedding)
% % == For more information about datasets (e.g. download links), please see the supplementary material 
% % randomLT_dataset.m : create dataset (documents with "random linear transformation" word embedding)
% % 
% % (b) Third-party toolboxes
% % figtreeKCenterClustering.m : the farthest-point clustering
% % BuildTreeMetric_HighDim_V2.m : clustering-based tree metric sampling
% % 
% % (c) Compute FlowAlign / DepthAlign function
% % L1_OT1D.m : univariate OT for weighted empirical measures (used for computing FlowAlign)
% %
% % --- A VIEW POINT OF EDGES OF PATHS FOR HIERARCHICAL ALIGNMENT IN TREE
% % DepthBasedTreeGWVE_L1_AlignedTree.m : compute DepthAlign 
% % == (with a view point on edges of paths for hierarchical alignment)
% % DepthBasedTreeGWVE_L1_AlignedTree_NR.m : compute DepthAlign
% % == (the nonrecursive version for computing DepthAlign)
% % == (with a view point on edges of paths for hierarchical alignment)
% %
% % --- A VIEW POINT OF VERTICES ON PATHS FOR HIERARCHICAL ALIGNMENT IN TREE 
% % Construct2DepthLevelTree: a preprocessing to construct 2-depth-level
% % trees for computing DepthAlign 
% % == based on a view point of vertices on paths from a root to each support for hierarchical alignment on tree 
% % L1_OT1D_Depth: compute OT for 2-depth-level tree (for "non-simple"
% % tree)
% % L1_OT1D_Depth_Simple: compute OT for 2-depth-level tree (where 1 tree
% % is "simple"
% % DepthBasedTreeGWVV_L1_AlignedTree.m : compute DepthAlign
% % == based on a view point of vertices on paths from a root to each support for hierarchical alignment on tree 
% % ------ For a variant: only considering mass in "leaves" for normalization
% % L1_OT1D_Depth_VX / L1_OT1D_Depth_Simple_VX / DepthBasedTreeGWVVX_L1_AlignedTree.m: are corresponding
% % functions of L1_OT1D_Depth / L1_OT1D_Depth_Simple / DepthBasedTreeGWVV_L1_AlignedTree.m
% %
% % OptimizeRoot_FA_L1OT: a version of FlowAlign where we optimize roots
% % using naive exhaustive search for root alignment.
% %
% % (d) Compute k-NN
% % compute_AlignedRoot_TSGW_Flow_DM2.m : compute pair-wise FlowAlign
% % compute_AlignedRoot_TSGW_DepthVE_DM2.m : compute pair-wise DepthAlign 
% % == based on a view point of edges on paths from a root to each support
% % for hierarchical alignment on tree.
% % compute_AlignedRoot_TSGW_DepthVV_DM2.m : compute pair-wise DepthAlign 
% % == based on a view point of vertices on paths from a root to each support
% % for hierarchical alignment on tree.
% % compute_AlignedRoot_TSGW_DepthVVX_DM2.m : compute pair-wise DepthAlign 
% % == based on a view point of vertices on paths from a root to each support
% % for hierarchical alignment on tree with a variant: only considering mass in "leaves" for normalization 
% % 
% % kNN2_TWITTER_AR_TSW_Flow.m : kNN with FlowAlign
% % kNN2_TWITTER_AR_TSW_DepthVE.m : kNN with DepthAlign
% % == a view point on edges of paths
% % kNN2_TWITTER_AR_TSW_DepthVV.m : kNN with DepthAlign
% % == a view point on vertices of paths
% % kNN2_TWITTER_AR_TSW_DepthVVX.m : kNN with DepthAlign
% % == a view point on vertices of paths with a variant for weight
% % normalization.
% % 
% % (*) Kmeans_FlowAlign_Barycenter
% % Code for K-means with FlowAlign (more details in the corresponding
% % folder)
% %
