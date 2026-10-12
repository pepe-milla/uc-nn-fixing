
% Reproducibility
rng(0);

% Load data: ND [n x 24], U (thermal, 16), Ug (turbine, 6), Up (pump, 6).
% Concatenate into a single commitment matrix: columns 1-16 thermal,
% 17-22 turbine, 23-28 pump.
load("../data/Training_Data.mat");
U = cat(2, U, Ug, Up);               % [n x 28 x 24]
numSamples = size(ND,1);
numGen  = size(U,2);                 % 28 units
numHours = size(U,3);                % 24 hours
idxTurb = 17:22;
idxPump = 23:28;
unitNames = [compose("G%d",(1:16)'); compose("T%d",(1:6)'); compose("B%d",(1:6)')];

% Groups of identical units: thermal pairs/trio as in the no-PHES case,
% plus the 6 identical turbine units and the 6 identical pump units.
% The MILP picks arbitrarily among identical units (same cost), so raw
% labels are degenerate; the stacked canonical encoding (unit j of a group
% is ON iff the group count is >= j) makes the demand->commitment relation
% a well-defined function, for both training and evaluation.
% Units 5 and 6 (the two CCGTs) are NOT identical (different Pmax, Mv and
% Pvc in the corrected data) and are therefore not grouped.
idGroups = {[1 2],[3 4],[7 8],[9 10],[11 12 13],[15 16], idxTurb, idxPump};
U = permute(canonicalize(permute(U,[3 2 1]), idGroups), [3 2 1]);

% Data partition: training (60%), validation (10%), calibration (15%) and
% test (15%). Validation is used ONLY for early stopping / checkpoint
% selection; since the model is selected to look good on it, its
% probabilities are optimistically biased and it must NOT be reused to
% calibrate thresholds or fixing bands. That is what the independent
% calibration set is for. Test remains untouched until the final metrics.
cv1 = cvpartition(numSamples, 'HoldOut', 0.40);
idxTrain = find(training(cv1));
idxHold  = find(test(cv1));
cv2 = cvpartition(numel(idxHold), 'HoldOut', 0.75);
idxVal   = idxHold(training(cv2));
idxRest  = idxHold(test(cv2));
cv3 = cvpartition(numel(idxRest), 'HoldOut', 0.5);
idxCal  = idxRest(training(cv3));
idxTest = idxRest(test(cv3));

prepX = @(idx) mat2cell(permute(ND(idx,:),[2 1]), numHours, ones(numel(idx),1));
prepY = @(idx) permute(U(idx,:,:), [3 2 1]);   % [24 x 28 x n]

XTrain = prepX(idxTrain);   YTrain = prepY(idxTrain);
XVal   = prepX(idxVal);     YVal   = prepY(idxVal);
XCal   = prepX(idxCal);     YCal   = prepY(idxCal);
XTest  = prepX(idxTest);    YTest  = prepY(idxTest);

% Neural network architecture (as described in the article): two stacked
% LSTM layers followed by fully connected layers. Note: a bidirectional
% LSTM (bilstmLayer) is arguably better suited to day-ahead UC (the
% commitment at hour t may depend on later hours of the known demand
% profile, even more so with the reservoir coupling all 24 hours); kept
% unidirectional for consistency with the article.
layers = [
    sequenceInputLayer(1,"Name","sequence","MinLength",numHours,"Normalization","zscore")
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

% Training options
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

% Training
tic
netTrained = trainnet(XTrain,YTrain,net,lossFcn,options);
t_train = toc;

% Per-unit decision thresholds calibrated on the CALIBRATION set
% (independent of training and of early stopping).
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

% Prediction on the test set
numTest = numel(XTest);
YPredBinary = false(numHours,numGen,numTest);
for i = 1:numTest
    YPred = predict(netTrained, XTest{i});
    YPredBinary(:,:,i) = YPred > thresholds;
end

% Assessment. Predictions are canonicalized too, so swapping two identical
% units does not count as an error (permutation-invariant metric).
YPredBinary = canonicalize(YPredBinary, idGroups);
Matches = YPredBinary == YTest;
accuracy1 = sum(Matches,"all")/numel(Matches);
accuracy2 = sum(all(Matches,2),[1 3])/(numTest*numHours);
accuracy3 = sum(all(Matches,1),[2 3])/(numTest*numGen);
accuracy4 = sum(all(Matches,[1 3]),2)/numGen;
accuracy5 = sum(all(Matches,[1 2]),3)/numTest;
fprintf("General accuracy: %2f %%\n", accuracy1*100);
fprintf("Hours with matches in all units: %2f %%\n", accuracy2*100)
fprintf("Units with matches in all hours of the day: %2f %%\n", accuracy3*100)
fprintf("Units with matches in all hours of all the days of the test sample: %2f %%\n", accuracy4*100)
fprintf("Days with matches in all units and hours: %2f %%\n", accuracy5*100);

% Per-unit accuracy (helps to see which units the net struggles with)
accPerGen = squeeze(mean(Matches,[1 3]));
fprintf("\nPer-unit accuracy (threshold):\n");
for g = 1:numGen
    fprintf("  %-3s: %6.2f %%  (th = %.2f)\n", unitNames(g), accPerGen(g)*100, thresholds(g));
end

% Physical-consistency diagnostic: a reversible unit cannot turbine and
% pump in the same hour, so predicted counts must satisfy nT + nP <= 6.
% (The canonical labels satisfy this by construction; thresholded
% predictions may not. Violations are handled at the fixing stage.)
nT = squeeze(sum(YPredBinary(:,idxTurb,:),2));
nP = squeeze(sum(YPredBinary(:,idxPump,:),2));
fprintf("\nHours with predicted turbine+pump count > 6: %.3f %%\n", ...
    mean(nT+nP > 6,"all")*100);

% ---- Local functions ----
function X = canonicalize(X, groups)
% Stacked (per-hour) canonical encoding for groups of identical units.
% X: [24 x numGen x n] (hours x units x days).
% Within a group, only the NUMBER of committed units per hour is determined
% by the problem; which unit covers it is arbitrary. Encoding: unit j of
% the group is ON at hour t iff the group count at hour t is >= j. Among
% all equivalent assignments, the stacked one minimizes start-ups.
for k = 1:numel(groups)
    idx = groups{k};
    cnt = sum(X(:,idx,:), 2);            % [24 x 1 x n] committed units
    for j = 1:numel(idx)
        X(:,idx(j),:) = cnt >= j;
    end
end
end
