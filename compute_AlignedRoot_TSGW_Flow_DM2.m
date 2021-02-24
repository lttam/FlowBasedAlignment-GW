%--------
% Compute FlowAlign
%--------

clear all
clc

nTS = 10; % number of slices

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
L1GW_Flow_SUM = zeros(N, N);
runTime_TM = 0;
runTime_dRX = 0;
runTime_L1GW_Flow = 0;

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

disp(['...compute TGW_Flow mapping']);
%%%%%%%%
XX_dRX = cell(N, 1);
WW_dRX = cell(N, 1);

% extension to general simplex vector (not only uniform vector) !!!
tic
for ii = 1:N
    node_XXII = XX_VertexID{ii};
    ww_XXII = WW{ii};
    
    dRX_XXII = zeros(length(node_XXII), 1);
    for jj = 1:length(dRX_XXII)
        % distance from root to each support
        dRX_XXII(jj) = sum(TM.Edge_Weight(TM.Vertex_EdgeIdPath{node_XXII(jj)}));
    end
    
    [XX_dRX{ii}, sorted_id_XXII] = sort(dRX_XXII);
    WW_dRX{ii} = ww_XXII(sorted_id_XXII);
end
runTime_dRX_ii = toc;
% accumulate
runTime_dRX = runTime_dRX + runTime_dRX_ii;

disp(['...compute L1GW_Flow']);
%%%%%%%
L1GW_Flow = zeros(N, N);
tic
for ii = 1:(N-1)
    if mod(ii, 50) == 0
        disp(['......' num2str(ii)]);
    end
    for jj=(ii+1):N
        % L1GW_Flow(ii, jj)        
        L1GW_Flow(ii, jj) = L1_OT1D(XX_dRX{ii}, XX_dRX{jj}, WW_dRX{ii}, WW_dRX{jj});
        L1GW_Flow(jj, ii) = L1GW_Flow(ii, jj);
    end
end
runTime_L1GW_Flow_ii = toc;
% accumulate
runTime_L1GW_Flow = runTime_L1GW_Flow + runTime_L1GW_Flow_ii;
L1GW_Flow_SUM = L1GW_Flow_SUM + L1GW_Flow;

end

% Average
L1GW_Flow = L1GW_Flow_SUM/nTS;

save(['twitter_RLT_AR_TSGW_Flow_DM2_K' num2str(TM_KC) 'L' num2str(TM_L) 'S' num2str(nTS) '.mat'], 'L1GW_Flow', ...
    'runTime_TM', 'runTime_dRX', 'runTime_L1GW_Flow', ...
    'TM_L', 'TM_KC', 'nTS');

disp('FINISH!');


