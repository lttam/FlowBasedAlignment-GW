function dd = L1_OT1D(XX, ZZ, wXX, wZZ)

% Computed a univariate OT for weighted empirical measures
% already sorted univariate supports

% l1 ground metric

dd = 0;

idXX = 1;
idZZ = 1;

constEPS = 1e-50;

while ((idXX <= length(XX)) && (idZZ <= length(ZZ)))
    if wXX(idXX) <= wZZ(idZZ)
        dd = dd + wXX(idXX)*abs(XX(idXX) - ZZ(idZZ));
        
        % update
        wZZ(idZZ) = wZZ(idZZ) - wXX(idXX);
        idXX = idXX + 1;
        if wZZ(idZZ) < constEPS
            idZZ = idZZ + 1;
        end
    else
        dd = dd + wZZ(idZZ)*abs(XX(idXX) - ZZ(idZZ));
        
        % update
        wXX(idXX) = wXX(idXX) - wZZ(idZZ);
        idZZ = idZZ + 1;
        if wXX(idXX) < constEPS
            idXX = idXX + 1;
        end
    end
end


