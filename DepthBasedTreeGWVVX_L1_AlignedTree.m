function dd = DepthBasedTreeGWVVX_L1_AlignedTree(TM, pathXX, pathZZ, weightXX, weightZZ, TTXX, TTZZ, idXX, idZZ, ww)

% (alignment hierarchically along vertices in paths from a root to each
% support -- based on a view point on vertices)
% a variant: only considering mass in "leaves" for normalization

% ---------------------------
% % % % % Input % % % % % % %
% ---------------------------

% TM: tree metric
% pathXX, weightXX: set of paths, each path is a set of
% edge (from root to each support) & their corresponding weight in weightXX
% --> same for pathZZ, weightZZ
% ==> using for computing the simple case !!!

% TTXX, TTZZ: constructed 2-depth-level tree from (set_path_XX, ww_path_XX,
% TM)

% idXX: the current node_id of XX in tree TM_X
% idZZ: the current node_id of ZZ in tree TM_Z

% ww: weight for this matching

% ---------------------------
% OUTPUT:
% ---------------------------
% dd: depth-based tree Gromov-Wasserstein variant

% FOR EXAMPLE
% start from the root --> id = 1, matching mass T = 1
% dd = DepthBasedTreeGW(..., TTXX, TTZZ, 1, 1, 1)

constEPS = 1e-20;
MAXQ = 50000;
queueProb = cell(MAXQ, 1);

dd = 0;

% QUEUE
prob.idXX = idXX;
prob.idZZ = idZZ;
prob.ww = ww;

queueProb{1} = prob;
length_queueProb = 1;

prob_ID = 1;

% % TTXX.vID = XX_vID(1:numNN);
% % TTXX.Tree_ID = XX_Tree_ID(1:numNN);
% % TTXX.Tree_DD = XX_Tree_DD(1:numNN);
% % TTXX.Tree_WW = XX_Tree_WW(1:numNN);
% % TTXX.SimpleFlag = XX_SimpleFlag(1:numNN);

while prob_ID <= length_queueProb

    % solve prob_ID
    ii_idXX = queueProb{prob_ID}.idXX;
    ii_idZZ = queueProb{prob_ID}.idZZ;
    ii_ww = queueProb{prob_ID}.ww;
    
    % find (DX, WX)
    ii_TTXX = find(TTXX.vID==ii_idXX); % ALWAYS FOUND!
    DX = TTXX.Tree_DD{ii_TTXX};
    WX = TTXX.Tree_WW{ii_TTXX};
    simpleFlagX = TTXX.SimpleFlag(ii_TTXX);
    setIDX = TTXX.Tree_ID{ii_TTXX};
    
    % find (DZ, WZ)
    ii_TTZZ = find(TTZZ.vID==ii_idZZ); % ALWAYS FOUND!
    DZ = TTZZ.Tree_DD{ii_TTZZ};
    WZ = TTZZ.Tree_WW{ii_TTZZ};
    simpleFlagZ = TTZZ.SimpleFlag(ii_TTZZ);
    setIDZ = TTZZ.Tree_ID{ii_TTZZ};
    
    % simpleFlagX + simpleFlagZ == 2 --> continue;
    
    % recursive case
    if (simpleFlagX + simpleFlagZ) == 0
        [ii_dd, ii_TT] = L1_OT1D_Depth_VX(DX, DZ, WX, WZ);
        dd = dd + ii_ww*ii_dd;
        ii_TT = ww*ii_TT;

        % Add new sub-problems into queue
        % only recursive for child nodes
        for ii_idx = 2:length(setIDX)
            prob.idXX = setIDX(ii_idx);
            for jj_idz = 2:length(setIDZ)
                prob.ww = ii_TT(ii_idx, jj_idz);
                if prob.ww >= constEPS
                    prob.idZZ = setIDZ(jj_idz);
                    % add new prob into queue
                    queueProb{length_queueProb + 1} = prob;
                    length_queueProb = length_queueProb + 1;
                end
            end
        end
    end
    
    %simple case
    if (simpleFlagX + simpleFlagZ) == 1
        if simpleFlagX == 1 % XX
            dd = dd + ww*L1_OT1D_Depth_Simple_VX(TM, pathZZ, weightZZ, ii_idZZ);
        else % ZZ
            dd = dd + ww*L1_OT1D_Depth_Simple_VX(TM, pathXX, weightXX, ii_idXX);
        end
    end
    
    prob_ID = prob_ID + 1;
end





