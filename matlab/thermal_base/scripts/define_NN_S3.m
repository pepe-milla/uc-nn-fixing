
% Strategy 3 (thermal case): the network input is augmented with the
% solution of the LP relaxation of the same UC instance, in addition to the
% 24-hour net demand profile. Everything else (architecture, canonical
% encoding, data partition, threshold calibration) is identical to
% define_NN.m, so that the comparison isolates the effect of the input.
%
% REQUIRED DATA: ../data/Training_Data_S3.mat with
%   ND [n x 24]        net demand
%   U  [n x 16 x 24]   commitment of the MIP solution (the label)
%   UR [n x 16 x 24]   commitment of the LP RELAXATION (in [0,1])

% Reproducibility
rng(0);

% ---- Load data ---------------------------------------------------------
load("../data/Training_Data_S3.mat", "ND", "U", "UR");
numSamples = size(ND,1);
numGen   = size(U,2);                % 16 thermal units
numHours = size(U,3);                % 24 hours

% Groups of identical units (units 5 and 6 are NOT identical: different
% Pmax, Mv and Pvc in the corrected data)
idGroups = {[1 2],[3 4],[7 8],[9 10],[11 12 13],[15 16]};

% Labels: stacked canonical encoding (unit j of a group is ON iff the group
% count is >= j), so permutations of identical units are not counted as
% different solutions.
U = permute(canonicalize(permute(U,[3 2 1]), idGroups), [3 2 1]);

% Relaxation: the same degeneracy affects the LP solution, so each group is
% sorted in descending order. Component j then means "the j-th largest
% relaxed commitment", the continuous counterpart of the stacked encoding.
UR = permute(sortGroups(permute(UR,[3 2 1]), idGroups), [3 2 1]);

% ---- Data partition: 60 / 10 / 15 / 15 ---------------------------------
cv1 = cvpartition(numSamples, 'HoldOut', 0.40);
idxTrain = find(training(cv1));
idxHold  = find(test(cv1));
cv2 = cvpartition(numel(idxHold), 'HoldOut', 0.75);
idxVal   = idxHold(training(cv2));
idxRest  = idxHold(test(cv2));
cv3 = cvpartition(numel(idxRest), 'HoldOut', 0.5);
idxCal  = idxRest(training(cv3));
idxTest = idxRest(test(cv3));

% Input normalization: the demand channel is z-scored with TRAINING-set
% statistics; the 16 relaxation channels are already in [0,1] and are left
% untouched. Done explicitly rather than with the input layer's built-in
% normalization because a relaxed channel that is constant over the
% training set (a base unit always at 1) would have zero variance.
ndMu    = mean(ND(idxTrain,:), "all");
ndSigma = std(ND(idxTrain,:), 0, "all");

prepX = @(idx) arrayfun(@(k) [ (ND(idx(k),:)'-ndMu)/ndSigma , ...
                               squeeze(UR(idx(k),:,:))' ], ...
                        (1:numel(idx))', 'UniformOutput', false);
prepY = @(idx) permute(U(idx,:,:), [3 2 1]);        % [24 x 16 x n]

XTrain = prepX(idxTrain);   YTrain = prepY(idxTrain);
XVal   = prepX(idxVal);     YVal   = prepY(idxVal);
XCal   = prepX(idxCal);     YCal   = prepY(idxCal);
XTest  = prepX(idxTest);    YTest  = prepY(idxTest);

numFeat = 1 + numGen;                               % 17 input channels

% ---- Network architecture ----------------------------------------------
% Same as define_NN.m (the article's network) except for the input width.
layers = [
    sequenceInputLayer(numFeat,"Name","sequence","MinLength",numHours,"Normalization","none")
    lstmLayer(128,"Name","lstm1","OutputMode","sequence")
    dropoutLayer(0.2,"Name","drop1")
    lstmLayer(64,"Name","lstm2","OutputMode","sequence")
    fullyConnectedLayer(64,"Name","fc1")
    reluLayer("Name","relu1")
    dropoutLayer(0.2,"Name","drop2")
    fullyConnectedLayer(32,"Name","fc2")
    reluLayer("Name","relu2")
    fullyConnectedLayer(numGen,"Name","fc_output")
    sigmoidLayer("Name","sigmoid")];
net = dlnetwork(layers);

% Loss function
lossFcn = "binary-crossentropy";

% ---- Training options --------------------------------------------------
miniBatchSize = 32;
valFreq = max(1, floor(numel(XTrain)/miniBatchSize));   % validate once per epoch
options = trainingOptions("adam", ...
    InitialLearnRate=2e-3, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropPeriod=25, ...
    LearnRateDropFactor=0.5, ...
    GradientThreshold=1, ...
    MiniBatchSize=miniBatchSize, ...
    MaxEpochs=150, ...
    Shuffle="every-epoch", ...
    Verbose=false, ...
    ValidationData={XVal, YVal}, ...
    ValidationFrequency=valFreq, ...
    ValidationPatience=10, ...
    OutputNetwork="best-validation-loss", ...
    Metrics="rmse", ...
    Plots="training-progress");

% ---- Training ----------------------------------------------------------
tic
netTrained = trainnet(XTrain,YTrain,net,lossFcn,options);
t_train = toc;

% ---- Per-unit decision thresholds (calibration set) --------------------
YCalProb = zeros(numHours,numGen,numel(XCal));
for i = 1:numel(XCal)
    YCalProb(:,:,i) = predict(netTrained, XCal{i});
end
thresholds = 0.5*ones(1,numGen);
for g = 1:numGen
    bestAcc = -inf;
    for th = 0.05:0.05:0.95
        acc = mean((squeeze(YCalProb(:,g,:)) > th) == squeeze(YCal(:,g,:)), "all");
        if acc > bestAcc
            bestAcc = acc;
            thresholds(g) = th;
        end
    end
end

% ---- Assessment on the test split --------------------------------------
numTest = numel(XTest);
YPredBinary = false(numHours,numGen,numTest);
for i = 1:numTest
    YPredBinary(:,:,i) = predict(netTrained, XTest{i}) > thresholds;
end
YPredBinary = canonicalize(YPredBinary, idGroups);
Matches = YPredBinary == YTest;
accuracy1 = sum(Matches,"all")/numel(Matches);
accuracy2 = sum(all(Matches,2),[1 3])/(numTest*numHours);
accuracy3 = sum(all(Matches,1),[2 3])/(numTest*numGen);
accuracy4 = sum(all(Matches,[1 3]),2)/numGen;
accuracy5 = sum(all(Matches,[1 2]),3)/numTest;
fprintf("General accuracy: %2f %%\n", accuracy1*100);
fprintf("Hours with matches in all generators: %2f %%\n", accuracy2*100)
fprintf("Generators with matches in all hours of the day: %2f %%\n", accuracy3*100)
fprintf("Generators with matches in all hours of all the days of the test sample: %2f %%\n", accuracy4*100)
fprintf("Days with matches in all generators and hours: %2f %%\n", accuracy5*100);

accPerGen = squeeze(mean(Matches,[1 3]));
fprintf("\nPer-generator accuracy (threshold):\n");
for g = 1:numGen
    fprintf("  G%-2d: %6.2f %%  (th = %.2f)\n", g, accPerGen(g)*100, thresholds(g));
end

% The normalization statistics must be reproduced at assessment time, so
% they travel with the network.
% save('../trained models/NN_S3.mat','netTrained','ndMu','ndSigma');

% ---- Local functions ---------------------------------------------------
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

function X = sortGroups(X, groups)
% Continuous counterpart of canonicalize: within each group of identical
% units the values are sorted in descending order at every hour, so that
% component j is the j-th largest. Makes the input invariant to the
% arbitrary labelling the LP relaxation gives to interchangeable units.
for k = 1:numel(groups)
    idx = groups{k};
    X(:,idx,:) = sort(X(:,idx,:), 2, 'descend');
end
end
