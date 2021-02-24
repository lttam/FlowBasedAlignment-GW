function [dd, pair, allDD] = OptimizeRoot_FA_L1OT(id1, id2, TM, XX_VertexID, XX_dRX_org, d11)

% optimize root for FlowAlign 
% using naive exhaustive search for optimal roots

dd = d11;
pair = [1; 1];

% compute all
all_X1_dRX_org = cell(TM.nVertices, 1);
all_X2_dRX_org = cell(TM.nVertices, 1);

X1_dRX = XX_dRX_org{id1}; % column vectors
X2_dRX = XX_dRX_org{id2};

% root (id=1)
all_X1_dRX_org{1} = X1_dRX;
all_X2_dRX_org{1} = X2_dRX;

% root (id=ii)
for ii = 2:TM.nVertices

    X1 = XX_VertexID{id1}; %row vectors
    X2 = XX_VertexID{id2};
    
    dRII = TreeMetricFromRoot(ii, TM);
    
    X1_dRXII = zeros(length(X1), 1); % new root --> supports (column)
    X2_dRXII = zeros(length(X2), 1);
    
    % X1
    for jj = 1:length(X1)
        X1_dRXII(jj) = dRII + X1_dRX(jj) - 2*TreeMetric_Root2CA(ii, X1(jj), TM);
    end
    
    % X2
    for jj = 1:length(X2)
        X2_dRXII(jj) = dRII + X2_dRX(jj) -2*TreeMetric_Root2CA(ii, X2(jj), TM);
    end
    
    all_X1_dRX_org{ii} = X1_dRXII;
    all_X2_dRX_org{ii} = X2_dRXII;
end

% sort
all_X1_dRX = cell(TM.nVertices, 1);
all_X2_dRX = cell(TM.nVertices, 1);

for ii = 1:TM.nVertices
    all_X1_dRX{ii} = sort(all_X1_dRX_org{ii});
    all_X2_dRX{ii} = sort(all_X2_dRX_org{ii});
end

% % dd = d11;
% % pair = [1; 1];
allDD = zeros(TM.nVertices, TM.nVertices);
for ii = 1:TM.nVertices
    wII = ones(length(all_X1_dRX{ii}), 1) / length(all_X1_dRX{ii});
    for jj = 1:TM.nVertices
        
        wJJ = ones(length(all_X2_dRX{jj}), 1) / length(all_X2_dRX{jj});
        allDD(ii, jj) = L1_OT1D(all_X1_dRX{ii}, all_X2_dRX{jj}, wII, wJJ);
        
        % compare
        if allDD(ii, jj) < dd
            dd = allDD(ii, jj);
            pair = [ii; jj];
        end
    end
end


end
