function TTXX = Construct2DepthLevelTree(XX, WW, TM)

% "preprocessing" for DepthAlign with vertex view
% (alignment hierarchically along vertices in paths from a root to each
% support -- based on a view point on vertices)

% % INPUT:

% CELL
% XX: set of paths --- each path is a set of edge_ids

% VECTOR
% WW: weight for each path

% % OUTPUT:
% TTXX.vID: array of root (vertices id) for 2-depth-level tree
% TTXX.Tree_DD: distance (root -- ID1, then concat with set_child_id)
% TTXX.Tree_WW: weight for each distance
% TTXX.Tree_ID: id of each node in tree
% TTXX.SimpleFlag: simple tree (1) otherwise (0)

% NOTE THAT:
% ID_Edge_HighNode = ID_EDGE + 1

% ROOT --> vID = 1;
MAXN = 20000;

XX_vID = zeros(MAXN, 1);
XX_Tree_DD = cell(MAXN, 1);
XX_Tree_WW = cell(MAXN, 1);
XX_Tree_ID = cell(MAXN, 1);
XX_SimpleFlag = zeros(MAXN, 1);

% FOR SAMPLED-TREE-METRIC TYPE
% --------root---------
numNN = 1;

vID = 1; % ROOT

set_ID = [vID, TM.Vertex_ChildId{vID}];
set_DD = zeros(length(set_ID), 1);
set_WW = zeros(length(set_ID), 1);

% child_id
for ii = 2:length(set_ID)
    child_id_ii = set_ID(ii);
    set_DD(ii) = TM.Edge_Weight(child_id_ii - 1);
    % for each path --> check edge_id = child_id_ii - 1
    for jj = 1:length(XX)
       if ismember(child_id_ii - 1, XX{jj})
            set_WW(ii) = set_WW(ii) + WW(jj);
       end
    end
end

XX_vID(numNN) = vID;
if sum(set_WW(2:end)) == 0
    XX_SimpleFlag(numNN) = 1;
    set_WW(1) = 1;
end

XX_Tree_ID{numNN} = set_ID;
XX_Tree_DD{numNN} = set_DD;
XX_Tree_WW{numNN} = set_WW;

% FOR EACH PATH
for iiPP = 1:length(XX)
    pathEdgeID_ii = XX{iiPP};
    
    % for each edge ID
    for jjPP = 1:length(pathEdgeID_ii)
        % node_id (deeper level) = edge_id + 1
        vID = pathEdgeID_ii(jjPP) + 1;
        if ismember(vID, XX_vID(1:numNN)) == 0
            % new node !!!
            numNN = numNN + 1;
            
            set_ID = [vID, TM.Vertex_ChildId{vID}];
            set_DD = zeros(length(set_ID), 1);
            set_WW = zeros(length(set_ID), 1);

            % child_id
            for ii = 2:length(set_ID)
                child_id_ii = set_ID(ii);
                set_DD(ii) = TM.Edge_Weight(child_id_ii - 1);
                % for each path --> check edge_id = child_id_ii - 1
                for jj = 1:length(XX)
                   if ismember(child_id_ii - 1, XX{jj})
                        set_WW(ii) = set_WW(ii) + WW(jj);
                   end
                   
                   % for root (vID)!!!
                   if ismember(vID - 1, XX{jj})
                        tmp = XX{jj};
                        if tmp(end) == (vID-1) % if leave
                            set_WW(1) = set_WW(1) + WW(jj);
                        end
                   end
                end
            end

            XX_vID(numNN) = vID;
            if sum(set_WW(2:end)) == 0
                XX_SimpleFlag(numNN) = 1;
            end

            XX_Tree_ID{numNN} = set_ID;
            XX_Tree_DD{numNN} = set_DD;
            XX_Tree_WW{numNN} = set_WW;
            
        end
    end
end

% return 1:numNN
TTXX.vID = XX_vID(1:numNN);
TTXX.Tree_ID = XX_Tree_ID(1:numNN);
TTXX.Tree_DD = XX_Tree_DD(1:numNN);
TTXX.Tree_WW = XX_Tree_WW(1:numNN);
TTXX.SimpleFlag = XX_SimpleFlag(1:numNN);

end