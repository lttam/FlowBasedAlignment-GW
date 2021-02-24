%--------
% K-NN for FlowAlign
%--------

clear all
clc

nTS = 10; % number of slices

TM_K = 4;
TM_L = 6;

load('twitter_RLT.mat');
% YY (prediction value)

load(['twitter_RLT_AR_TSGW_Flow_DM2_K' num2str(TM_K) 'L' num2str(TM_L) 'S' num2str(nTS) '.mat']);
% L1GW_Flow

perTrain = 0.8;
NTr = round(N*perTrain);
NTe = N - NTr;

n_repeat = 20;
maxK = 30;

arrayACC = zeros(n_repeat, maxK);

for ii_repeat = 1:n_repeat

    disp(['......' num2str(ii_repeat)]);
    
    randID = randperm(N);
    IDTr = randID(1:NTr);
    IDTe = randID((NTr+1):end);
    
    YYTr = YY(IDTr);
    YYTe = YY(IDTe);

    XXL1 = L1GW_Flow(IDTe, IDTr);
    
    % prediction
    preYL1 = zeros(length(IDTe), maxK);
    
    for ii = 1:length(IDTe)
        
        % L1 
        [valL1, idL1] = sort(XXL1(ii, :));
        
        % sorted YY
        YYL1 = YYTr(idL1);
        
        % average
        preYL1(ii, 1) = YYL1(1);
        for kk = 2:maxK        
            preYL1(ii, kk) = mode(YYL1(1:kk));
        end
    end
    
    % MAE
    for kk = 1:maxK
        % L1 -- column vector
        arrayACC(ii_repeat, kk) = sum(abs(YYTe - preYL1(:, kk))==0) / NTe; 
    end
end

meanACC = mean(arrayACC);
stdACC = std(arrayACC);

save(['twitter_RLT_AR_TSGW_Flow_Results2_K' num2str(TM_K) 'L' num2str(TM_L) 'S' num2str(nTS) '.mat'], ...
    'arrayACC', 'meanACC', 'stdACC', ...
    'perTrain', 'n_repeat', 'maxK');

disp('FINISH!');









