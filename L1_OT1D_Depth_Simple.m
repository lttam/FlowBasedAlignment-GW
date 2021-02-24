function dd = L1_OT1D_Depth_Simple(TM, pathCell, WW, id)

LL = zeros(length(WW), 1);

% for each path
for ii = 1:length(pathCell)
    pathEdgeID = pathCell{ii};
    
    % ID_EDGE = ID_NODE(HIGH) - 1
    id_II = find(pathEdgeID == (id-1));
    if isempty(id_II) == 0 % found at id_II
        edgeLengthPath = TM.Edge_Weight(pathEdgeID);
        % length of path from id --> leaves
        LL(ii) = sum(edgeLengthPath((id_II+1):end)); 
    end
end

idNZ = find(LL~=0);

LL_NZ = LL(idNZ);

WW_NZ = WW(idNZ);
WW_NZ = WW_NZ / sum(WW_NZ); % normalization for weights

dd = sum(LL_NZ .* WW_NZ);

end


