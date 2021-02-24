function dd = SL2_OT1D(XX, ZZ, wXX, wZZ)

% Computed a univariate OT for weighted empirical measures
% already sorted univariate supports

% squared L2 ground metric

dd = 0;

idXX = 1;
idZZ = 1;

constEPS = 1e-50;

while ((idXX <= length(XX)) && (idZZ <= length(ZZ)))
    if wXX(idXX) <= wZZ(idZZ)
        dd = dd + wXX(idXX)*(XX(idXX) - ZZ(idZZ))^2;
        
        % update
        wZZ(idZZ) = wZZ(idZZ) - wXX(idXX);
        idXX = idXX + 1;
        if wZZ(idZZ) < constEPS
            idZZ = idZZ + 1;
        end
    else
        dd = dd + wZZ(idZZ)*(XX(idXX) - ZZ(idZZ))^2;
        
        % update
        wXX(idXX) = wXX(idXX) - wZZ(idZZ);
        idZZ = idZZ + 1;
        if wXX(idXX) < constEPS
            idXX = idXX + 1;
        end
    end
end


