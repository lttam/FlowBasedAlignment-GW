function [dd, orgTT] = L1_OT1D_Depth_VX(DX, DZ, WX, WZ)

% a variant: only considering mass in "leaves" for normalization

dd = 0;

% normalize weight
% dismiss root at id1
WX = WX / sum(WX(2:end));
WZ = WZ / sum(WZ(2:end));

% sort
% not affected
[XX, sortDX_ID] = sort(DX);
[ZZ, sortDZ_ID] = sort(DZ);

wXX = WX(sortDX_ID);
wZZ = WZ(sortDZ_ID);

% dismiss root at id1
idXX = 2;
idZZ = 2;

TT = zeros(length(wXX), length(wZZ));
constEPS = 1e-50;

while ((idXX <= length(XX)) && (idZZ <= length(ZZ)))
    if wXX(idXX) <= wZZ(idZZ)
        dd = dd + wXX(idXX)*abs(XX(idXX) - ZZ(idZZ));
        TT(idXX, idZZ) = wXX(idXX);
        
        % update
        wZZ(idZZ) = wZZ(idZZ) - wXX(idXX);
        idXX = idXX + 1;
        if wZZ(idZZ) < constEPS
            idZZ = idZZ + 1;
        end
    else
        dd = dd + wZZ(idZZ)*abs(XX(idXX) - ZZ(idZZ));
        TT(idXX, idZZ) = wZZ(idZZ);
        
        % update
        wXX(idXX) = wXX(idXX) - wZZ(idZZ);
        idZZ = idZZ + 1;
        if wXX(idXX) < constEPS
            idXX = idXX + 1;
        end
    end
end

% return the original index
[~, orgIDX] = sort(sortDX_ID);
[~, orgIDZ] = sort(sortDZ_ID);

orgTT = TT(orgIDX, orgIDZ); 




