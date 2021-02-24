%--------
% Compute DepthAlign
% (alignment hierarchically along vertices in paths from a root to each
% support -- based on a view point on vertices)
%--------

clear all
clc

nTS = 1; % number of slices

% build tree metric for the qm7 dataset
% molecules are centered to (0, 0, 0)
% build tree metric (same structure) for each molecules from
%       all data positions
% compute flow-based representation (support --> distance from
%       root to each support)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% using clustering-based tree metric construction
% from tree-sliced TW

load('twitter_RLT.mat');

% function [TM, XX_VertexID] = BuildTreeMetric_HighDim_V2(XX, L, KC)
%%%%%%%%
% INPUT:
% XX: cell of supports for N empirical measures
% Each element XX{ii} is a matrix of supports (N_i x dim) where N_i is the
% number of supports with dim dimensions.
% L: the predefined highest level of tree
% (as a stopping condition)
% KC: number of clusters (for the farthest clustering of Gonzalez)
% (when one partitions space by clustering)
%%%%%%%%
% OUTPUT:
% TM: Tree structure T
% XX_VertexID: corresponding for XX in tree.

TM_L = 6; %highest level
TM_KC = 4; %# of clusters

% OUTPUT
L1GW_Depth_SUM = zeros(N, N);
runTime_TM = 0;
runTime_Mapping_Depth = 0;
runTime_L1GW_Depth = 0;

% REPEAT for each tree slice
for idSS = 1:nTS

    disp(['.........Tree: #' num2str(idSS)]);
    
%%%%%%%%
disp(['...compute tree metric']);
tic
[TM, XX_VertexID] = BuildTreeMetric_HighDim_V2(XX_center, TM_L, TM_KC);
runTime_TM_ii = toc;
% accumulate
runTime_TM = runTime_TM + runTime_TM_ii;

disp(['...compute TGW_Depth mapping']);
%%%%%%%%
XX_path_cell = cell(N, 1);
% % XX_weight_cell = cell(N, 1);
XX_weight_cell = WW;
tic
for ii = 1:N
    node_XXII = XX_VertexID{ii};
    
    path_XXII = cell(length(node_XXII), 1);
    for jj = 1:length(path_XXII)
        % distance from root to each support
        path_XXII{jj} = TM.Vertex_EdgeIdPath{node_XXII(jj)};
    end
    
    XX_path_cell{ii} = path_XXII;
    
% %     % WEIGHT FOR EACH PATH --> OR EACH SUPPORT (ALREADY !!!)
% %     XX_weight_cell{ii} = 1.0*ones(length(path_XXII), 1)/length(path_XXII);
 
end
runTime_Mapping_Depth_ii = toc;
% accumulate
runTime_Mapping_Depth = runTime_Mapping_Depth + runTime_Mapping_Depth_ii;

disp(['......compute TGW_Depth mapping --> construct 2-depth-level trees']);
TTXX = cell(N, 1);
for ii = 1:N
    TTXX{ii} = Construct2DepthLevelTree(XX_path_cell{ii}, XX_weight_cell{ii}, TM);
end

disp(['...compute L1GW_Depth']);
%%%%%%%
L1GW_Depth = zeros(N, N);
tic
for ii = 1:(N-1)
    if mod(ii, 10) == 0
        disp(['......' num2str(ii)]);
    end
    for jj=(ii+1):N

        L1GW_Depth(ii, jj) = DepthBasedTreeGWVV_L1_AlignedTree(TM, XX_path_cell{ii}, XX_path_cell{jj}, ...
            XX_weight_cell{ii}, XX_weight_cell{jj}, TTXX{ii}, TTXX{jj}, 1, 1, 1);

        L1GW_Depth(jj, ii) = L1GW_Depth(ii, jj);
    end
end
runTime_L1GW_Depth_ii = toc;
% accumulate
runTime_L1GW_Depth = runTime_L1GW_Depth + runTime_L1GW_Depth_ii;
L1GW_Depth_SUM = L1GW_Depth_SUM + L1GW_Depth;

end

% Average
L1GW_Depth = L1GW_Depth_SUM/nTS;

save(['twitter_RLT_AR_TSGW_DepthVV_DM2_K' num2str(TM_KC) 'L' num2str(TM_L) 'S' num2str(nTS) '.mat'], 'L1GW_Depth', ...
    'runTime_TM', 'runTime_L1GW_Depth', 'runTime_Mapping_Depth', ...
    'TM_L', 'TM_KC', 'nTS');

disp('FINISH!');




