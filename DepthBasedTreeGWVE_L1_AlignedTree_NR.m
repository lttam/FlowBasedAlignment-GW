function dd = DepthBasedTreeGWVE_L1_AlignedTree_NR(XX, ZZ, TM, pXX, pZZ)

% NON recursive version!!!

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
% ---------------------------
% OUTPUT:
% ---------------------------
% dd: depth-based tree Gromov-Wasserstein variant

% FOR EXAMPLE
% dd = DepthBasedTreeGW(XX, ZZ, TMX, TMZ, pXX, pZZ, 1)


% % XX, ZZ, TM, pXX, pZZ, ww
% initialization with the given pair in the stack
allXX{1} = XX;
allZZ{1} = ZZ;
allPX{1} = pXX;
allPZ{1} = pZZ;
allWW(1) = 1;
allID(1) = 1;

dd = 0;
while ~isempty(allXX)
    % compute for the first pair
    % collect the first pair
    XX = allXX{1};
    ZZ = allZZ{1};
    pXX = allPX{1};
    pZZ = allPZ{1};
    ww = allWW(1);
    id = allID(1);
    
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

    % ignore 0-0 case
    if sum_length_XX == 0
        % solve 0-X case (simple case)
        % check ZZ!
        if sum_length_ZZ > 0
            wZZ = ww / sum(pZZ); % weight for each support in ZZ;
            for ii = 1:length(ZZ) % each support path (edge id)
                % path from the (id)^th level --> support
                ddZI = sum(TM.Edge_Weight(ZZ{ii}));
                % weighted by local pZZ(ii)
                % then weighted by global wZZ
                ddZI = wZZ*ddZI*pZZ(ii);
                dd = dd + ddZI;
            end
        end
    elseif sum_length_ZZ == 0
        % solve X-0 case (simple case)
        % check XX!
        if sum_length_XX > 0
            wXX = ww / sum(pXX); % weight for each support in XX;
            for ii = 1:length(XX) % each support path (edge id)
                % path from the (id)^th level --> support
                ddXI = sum(TM.Edge_Weight(XX{ii}));
                % weighted by local pXX(ii)
                % then weighted by global wXX
                ddXI = wXX*ddXI*pXX(ii);
                dd = dd + ddXI;
            end
        end
    else
        % solve X-X case (recursive case without recursive implementation)
        % solve the 2-depth-level tree & put recursive cases into the queue
        
        % compute distance 2-depth-level tree !!!
        % at the (id)^th order
        XX_ID = zeros(length(XX), 1);
        % extract elements having length equal or greater than id
        iiGID = find(cellfun(@(c) length(c), XX) >= id);
        % update value for XX_ID
        XX_ID(iiGID) = cellfun(@(c) TM.Edge_Weight(c(id)), XX(iiGID));
        
        % same for ZZ
        ZZ_ID = zeros(length(ZZ), 1);
        % extract elements having length equal or greater than id
        iiGID = find(cellfun(@(c) length(c), ZZ) >= id);
        % update value for XX_ID
        ZZ_ID(iiGID) = cellfun(@(c) TM.Edge_Weight(c(id)), ZZ(iiGID));
        
        % global weights
        wXX = ww / sum(pXX);
        wZZ = ww / sum(pZZ);

        % sort & diff --> sorted unique (row)
        sorted_supp_XX = sort(XX_ID);
        sorted_supp_XX = sorted_supp_XX([true; diff(sorted_supp_XX(:))>0]);
        % for ZZ_ID
        sorted_supp_ZZ = sort(ZZ_ID);
        sorted_supp_ZZ = sorted_supp_ZZ([true; diff(sorted_supp_ZZ(:))>0]);
        
        % column vector
        sorted_supp_XX = sorted_supp_XX(:);
        sorted_supp_ZZ = sorted_supp_ZZ(:);

        % corresponding sorted weight (following the order of supports)
        sorted_ww_XX = zeros(length(sorted_supp_XX), 1);
        path_XX_cell = cell(length(sorted_ww_XX), 1);
        weight_XX_cell = cell(length(sorted_ww_XX), 1);
        for ii = 1:length(sorted_supp_XX)
            idII = find(XX_ID == sorted_supp_XX(ii));
            % save
            path_XX_cell{ii} = XX(idII);
            weight_XX_cell{ii} = pXX(idII);

            sorted_ww_XX(ii) = wXX*sum(pXX(idII));
        end

        % same for ZZ
        sorted_ww_ZZ = zeros(length(sorted_supp_ZZ), 1);
        path_ZZ_cell = cell(length(sorted_ww_ZZ), 1);
        weight_ZZ_cell = cell(length(sorted_ww_ZZ), 1);
        for ii = 1:length(sorted_supp_ZZ)
            idII = find(ZZ_ID == sorted_supp_ZZ(ii));
            % save
            path_ZZ_cell{ii} = ZZ(idII);
            weight_ZZ_cell{ii} = pZZ(idII);

            sorted_ww_ZZ(ii) = wZZ*sum(pZZ(idII));
        end
     
        % compute 1D-OT (sorted_supp_XX, sorted_ww_XX) or its corresponding
        % with ZZ

        % recursive info:
        % %     path_ZZ_cell 
        % %     weight_ZZ_cell 

        ii_id_XX = 1;
        ii_id_ZZ = 1;
        while ((ii_id_XX <= length(sorted_ww_XX)) && (ii_id_ZZ <= length(sorted_ww_ZZ)))

            if sorted_ww_XX(ii_id_XX) <= sorted_ww_ZZ(ii_id_ZZ)

                dd = dd + sorted_ww_XX(ii_id_XX)*abs(sorted_supp_XX(ii_id_XX) - sorted_supp_ZZ(ii_id_ZZ));
                
                % add new pairs to queue!!!
                allXX{end+1} = path_XX_cell{ii_id_XX};
                allZZ{end+1} = path_ZZ_cell{ii_id_ZZ};
                allPX{end+1} = weight_XX_cell{ii_id_XX};
                allPZ{end+1} = weight_ZZ_cell{ii_id_ZZ};
                allWW(end+1) = sorted_ww_XX(ii_id_XX);
                allID(end+1) = id+1;

                % update
                sorted_ww_ZZ(ii_id_ZZ) = sorted_ww_ZZ(ii_id_ZZ) - sorted_ww_XX(ii_id_XX);
                ii_id_XX = ii_id_XX + 1;
                if sorted_ww_ZZ(ii_id_ZZ) < 1e-50 % compare to zero
                    ii_id_ZZ = ii_id_ZZ + 1;
                end

            else

                dd = dd + sorted_ww_ZZ(ii_id_ZZ)*abs(sorted_supp_XX(ii_id_XX) - sorted_supp_ZZ(ii_id_ZZ));
 
                % add new pairs to queue!!!
                allXX{end+1} = path_XX_cell{ii_id_XX};
                allZZ{end+1} = path_ZZ_cell{ii_id_ZZ};
                allPX{end+1} = weight_XX_cell{ii_id_XX};
                allPZ{end+1} = weight_ZZ_cell{ii_id_ZZ};
                allWW(end+1) = sorted_ww_ZZ(ii_id_ZZ);
                allID(end+1) = id+1;
    
                % update
                sorted_ww_XX(ii_id_XX) = sorted_ww_XX(ii_id_XX) - sorted_ww_ZZ(ii_id_ZZ);
                ii_id_ZZ = ii_id_ZZ + 1;
                if sorted_ww_XX(ii_id_XX) < 1e-50 % compare to zero
                    ii_id_XX = ii_id_XX + 1;
                end
            end
        end
    end

    
    % delete the first pair (already computed)
    allXX(1) = [];
    allZZ(1) = [];
    allPX(1) = [];
    allPZ(1) = [];
    allWW(1) = [];
    allID(1) = [];
end

end
       
     


