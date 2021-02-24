%--------
% generate twitter dataset with "random linear transformation" for word
% embeddings
%--------

clear all
clc

% load original dataset (preprocessed document datasets with word
% embeddings)
load('twitter_org.mat');

N = length(X);

XX = cell(N, 1); % each cell: (N x dim) matrix
WW = cell(N, 1); % each cell: column vector

XX_center = cell(N, 1); % each cell: centerized by 0-vector
XX_DM = cell(N, 1); % each cell: L2-distance matrix

YY = Y';

% randomly linear transformation for word-embedding
dim = size(X{1}, 1);

for ii = 1:N
    
    % normalization for weights
    WW{ii} = BOW_X{ii}'/sum(BOW_X{ii});
    
    % randomly linear transformation matrix
    RM = randn(dim, dim);
    XI_RM = RM * X{ii}; %dim x N
    XI_center = XI_RM - repmat(mean(XI_RM, 2), 1, size(XI_RM, 2)); %dim x N
    
    XX{ii} = XI_RM';
    XX_center{ii} = XI_center';
     
end

save('twitter_RLT.mat', 'N', 'XX', 'XX_center', 'YY', 'WW');

disp('FINISH!!!');


