% Strategy 1 (thermal case): identification of the "trustful units".
% Standalone script: it only needs the saved network and the training-period
% dataset, so it can be run at any time after training.
%
% Criterion (as in the article): a unit is trustful if its on/off pattern is
% predicted without a single error at every hour of every day of the
% TRAINING PERIOD (2019-2021). No day of that period belongs to the test
% period, so the selection never looks at the test set.
%
% Note on canonical encoding: after the stacked encoding a "unit" is a layer
% of its group ("at least j units of the group are ON"), so a trustful layer
% means that level of the group count is perfectly predictable.

dataFile = '../data/Training_Data_2019-2021.mat';
netFile  = '../trained models/NN.mat';
outFile  = '../data/Strategy1_units.mat';

% ---- Load data and network ---------------------------------------------
load(dataFile, 'ND', 'U');
load(netFile, 'netTrained');

numSamples = size(ND,1);
numGen     = size(U,2);              % 16 thermal units
numHours   = size(U,3);
unitNames  = compose("G%d",(1:numGen)');

% Same groups as in define_NN.m (units 5 and 6 are NOT identical)
idGroups = {[1 2],[3 4],[7 8],[9 10],[11 12 13],[15 16]};

% ---- Targets and predictions over the whole training period ------------
Y = canonicalize(permute(U, [3 2 1]), idGroups);        % [24 x numGen x n]

YPred = false(numHours, numGen, numSamples);
for d = 1:numSamples
    YPred(:,:,d) = nnPredict(netTrained, ND(d,:)', numGen) > 0.5;
end
YPred = canonicalize(YPred, idGroups);

accTrain = squeeze(mean(YPred == Y, [1 3]));
trustful = accTrain(:)' == 1;
onShare  = squeeze(mean(Y, [1 3]));

% ---- Report ------------------------------------------------------------
fprintf("Strategy 1: trustful units (100%% accuracy over the %d days of the training period)\n\n", numSamples);
fprintf("Unit | accuracy | trustful | ON-share\n");
fprintf("-----+----------+----------+---------\n");
for g = 1:numGen
    if trustful(g), mark = "  yes   "; else, mark = "   no   "; end
    fprintf("%-4s |  %6.4f  | %s |  %5.1f %%\n", ...
        unitNames(g), accTrain(g), mark, onShare(g)*100);
end
fprintf("\nTrustful units: %d of %d (%.1f %%) -> %d of %d binaries/day\n", ...
    nnz(trustful), numGen, 100*nnz(trustful)/numGen, ...
    nnz(trustful)*numHours, numGen*numHours);

trivial = trustful & (onShare(:)' == 0 | onShare(:)' == 1);
if any(trivial)
    fprintf("Of these, %d are constant throughout the training period (%s)\n", ...
        nnz(trivial), strjoin(cellstr(unitNames(trivial)), ", "));
end

save(outFile, 'trustful', 'accTrain', 'onShare', 'idGroups', 'unitNames');
fprintf("\nSaved trustful-unit mask to %s\n", outFile);

% ---- Local functions ---------------------------------------------------
function P = nnPredict(net, x, numGen)
% Inference compatible with both older (formatted-dlarray only) and newer
% MATLAB releases. Returns [24 x numGen].
y = extractdata(predict(net, dlarray(single(x), 'TC')));
if size(y,1) == numGen, y = y'; end
P = double(y);
end

function X = canonicalize(X, groups)
% Stacked (per-hour) canonical encoding for groups of identical units:
% unit j of a group is ON at hour t iff the group count at hour t is >= j.
for k = 1:numel(groups)
    idx = groups{k};
    cnt = sum(X(:,idx,:), 2);
    for j = 1:numel(idx)
        X(:,idx(j),:) = cnt >= j;
    end
end
end
