function dd = DepthBasedTreeGWVE_L1_AlignedTree(XX, ZZ, TM, pXX, pZZ, ww, id)

% (alignment hierarchically over edges along paths from a root to each
% support -- based on a view point on edges)

% ---------------------------
% % % % % Input % % % % % % %
% ---------------------------
% XX: cell of all paths from root to each support (edge-id path) in tree
% metric TMX

% XX: cell of all paths from root to each support (edge-id path) in tree
% metric TMX

% TMXX, TMZZ: tree metrics

% ww: weight for this matching
% id: the considered order of the path

% ---------------------------
% OUTPUT:
% ---------------------------
% dd: depth-based tree Gromov-Wasserstein variant

% FOR EXAMPLE
% dd = DepthBasedTreeGWVE_L1_AlignedTree(XX, ZZ, TMX, TMZ, pXX, pZZ, 1, 1)

if ww <= 0
    dd = 0;
else

% check empty of XX
sum_length_XX = 0;
for ii = 1:length(XX)
    % length of XX{ii} at the (id)^th order
    sum_length_XX = sum_length_XX + (length(XX{ii}) - id + 1);
    if sum_length_XX > 0
        break;
    end
end

% check empty of ZZ
sum_length_ZZ = 0;
for ii = 1:length(ZZ)
    % length of XX{ii} at the (id)^th order
    sum_length_ZZ = sum_length_ZZ + (length(ZZ{ii}) - id + 1);
    if sum_length_ZZ > 0
        break;
    end
end

if sum_length_XX == 0
    
    % compute moving for ZZ
    dd = 0;
    if sum_length_ZZ > 0
        wZZ = ww / sum(pZZ); % weight for each support in ZZ;
        for ii = 1:length(ZZ) % each support path (edge id)
            % path from the (id)^th level --> support
            ZI = ZZ{ii}; % edge-id path
            ddZI = 0;
            for jj = id:length(ZI)
                ddZI = ddZI + TM.Edge_Weight(ZI(jj));
            end
            ddZI = wZZ*ddZI*pZZ(ii);
            dd = dd + ddZI;
        end
    end
    
elseif sum_length_ZZ == 0
    
    % compute moving for XX
    dd = 0;
    if sum_length_XX > 0
        wXX = ww / sum(pXX); % weight for each support in ZZ;
        for ii = 1:length(XX) % each support path (edge id)
            % path from the (id)^th level --> support
            XI = XX{ii}; % edge-id path
            ddXI = 0;
            for jj = id:length(XI)
                ddXI = ddXI + TM.Edge_Weight(XI(jj));
            end
            ddXI = wXX*ddXI*pXX(ii);
            dd = dd + ddXI;
        end
    end
    
else
    
    % solve & call recursive steps
    
    % compute distance 2-depth-level tree !!!
    % at the (id)^th order
    XX_ID = zeros(length(XX), 1);
    for ii = 1:length(XX)
        XI = XX{ii};
        if length(XI) < id
            XX_ID(ii) = 0;
        else
            XX_ID(ii) = TM.Edge_Weight(XI(id));
        end
    end
    % same for ZZ
    ZZ_ID = zeros(length(ZZ), 1);
    for ii = 1:length(ZZ)
        ZI = ZZ{ii};
        if length(ZI) < id
            ZZ_ID(ii) = 0;
        else
            ZZ_ID(ii) = TM.Edge_Weight(ZI(id));
        end
    end
    
    wXX = ww / sum(pXX);
    wZZ = ww / sum(pZZ);
    
    sorted_supp_XX = sort(unique(XX_ID));
    sorted_supp_ZZ = sort(unique(ZZ_ID));
    % column vector
    sorted_supp_XX = sorted_supp_XX(:);
    sorted_supp_ZZ = sorted_supp_ZZ(:);
    
    % corresponding sorted weight (following the order of supports)
    sorted_ww_XX = zeros(length(sorted_supp_XX), 1);
    path_XX_cell = cell(length(sorted_ww_XX), 1);
    weight_XX_cell = cell(length(sorted_ww_XX), 1);
    for ii = 1:length(sorted_supp_XX)
        idII = find(XX_ID == sorted_supp_XX(ii));
        
        path_XXII = XX(idII); % cell
        weight_XXII = pXX(idII); % vector
        % save
        path_XX_cell{ii} = path_XXII;
        weight_XX_cell{ii} = weight_XXII;
        
        sorted_ww_XX(ii) = wXX*sum(weight_XXII);
    end
    
    % same for ZZ
    sorted_ww_ZZ = zeros(length(sorted_supp_ZZ), 1);
    path_ZZ_cell = cell(length(sorted_ww_ZZ), 1);
    weight_ZZ_cell = cell(length(sorted_ww_ZZ), 1);
    for ii = 1:length(sorted_supp_ZZ)
        idII = find(ZZ_ID == sorted_supp_ZZ(ii));
        
        path_ZZII = ZZ(idII); % cell
        weight_ZZII = pZZ(idII); % vector
        % save
        path_ZZ_cell{ii} = path_ZZII;
        weight_ZZ_cell{ii} = weight_ZZII;
        
        sorted_ww_ZZ(ii) = wZZ*sum(weight_ZZII);
    end
    
    % compute 1D-OT (sorted_supp_XX, sorted_ww_XX) or its corresponding
    % with ZZ
    
    % recursive info:
    % %     path_ZZ_cell 
    % %     weight_ZZ_cell 
    
    ii_id_XX = 1;
    ii_id_ZZ = 1;
    dd = 0;
    while ((ii_id_XX <= length(sorted_ww_XX)) && (ii_id_ZZ <= length(sorted_ww_ZZ)))
    
        if sorted_ww_XX(ii_id_XX) <= sorted_ww_ZZ(ii_id_ZZ)
            
            dd = dd + sorted_ww_XX(ii_id_XX)*abs(sorted_supp_XX(ii_id_XX) - sorted_supp_ZZ(ii_id_ZZ));
            % recursive
            dd = dd + DepthBasedTreeGWVE_L1_AlignedTree(path_XX_cell{ii_id_XX}, path_ZZ_cell{ii_id_ZZ}, ...
                TM, weight_XX_cell{ii_id_XX}, weight_ZZ_cell{ii_id_ZZ}, sorted_ww_XX(ii_id_XX), id + 1);
            
            % update
            sorted_ww_ZZ(ii_id_ZZ) = sorted_ww_ZZ(ii_id_ZZ) - sorted_ww_XX(ii_id_XX);
            ii_id_XX = ii_id_XX + 1;
            if sorted_ww_ZZ(ii_id_ZZ) < 1e-50 % compare to zero
                ii_id_ZZ = ii_id_ZZ + 1;
            end
            
        else
            
            dd = dd + sorted_ww_ZZ(ii_id_ZZ)*abs(sorted_supp_XX(ii_id_XX) - sorted_supp_ZZ(ii_id_ZZ));
            % recursive
            dd = dd + DepthBasedTreeGWVE_L1_AlignedTree(path_XX_cell{ii_id_XX}, path_ZZ_cell{ii_id_ZZ}, ...
                TM, weight_XX_cell{ii_id_XX}, weight_ZZ_cell{ii_id_ZZ}, sorted_ww_ZZ(ii_id_ZZ), id + 1);
            
            % update
            sorted_ww_XX(ii_id_XX) = sorted_ww_XX(ii_id_XX) - sorted_ww_ZZ(ii_id_ZZ);
            ii_id_ZZ = ii_id_ZZ + 1;
            if sorted_ww_XX(ii_id_XX) < 1e-50 % compare to zero
                ii_id_XX = ii_id_XX + 1;
            end
            
        end
        
    end
    
end

end
    






