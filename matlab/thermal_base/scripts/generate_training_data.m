
% Initial and final day of the training period
init_time = datetime('2019-12-19 20:00','InputFormat','yyyy-MM-dd HH:mm');
final_time = datetime('2021-12-19 20:00','InputFormat','yyyy-MM-dd HH:mm');

% Baseline guard: the model includes data_gms/fixed_variables.gms
% unconditionally, so any content left there by a previous NN campaign
% would silently constrain every instance solved here. Emptied once,
% before the loop.
modelDir = fullfile(uc_config().gamsRoot, 'thermal_base');
fid = fopen(fullfile(modelDir,'data_gms','fixed_variables.gms'), 'w');
fprintf(fid, '* No fixed variables (label generation)\n');
fclose(fid);
[DAYDEM, DAYWIND] = uc_load_profiles('2019-2021');
numDays = days(duration(final_time-init_time));
DayDemand = zeros(1,24);
DayWind = zeros(1,24);

% Result variables for storage
Z = zeros(numDays,1);
U = zeros(numDays,16,24);
DELTA = zeros(numDays,16,24);
POWER = zeros(numDays,16,24);
WC = zeros(numDays,24);
COMP_TIME = zeros(numDays,1);
OPT_GAP = zeros(numDays,1);
T = zeros(numDays,1);
ND = zeros(numDays,24);


for d = 1:numDays
    % Perform the average of the 5min-data in order to get the hourly
    % demand and wind power production.
    DayDemand = DAYDEM(d,:);
    DayWind   = DAYWIND(d,:);
    % Store net demand
    ND(d,:) = DayDemand - DayWind;

    % Write the demand and wind file in GAMS environment 
    write_net_demand;

    % Run GAMS
    cd('../')
    run_file='./MIP_GCI_SIMPL_sinPHSP.gms';
    ejecuta=['"D:\Program Files\GAMS\GAMS23.5\gams.exe"' ' ' run_file ' lo=2'];
    % Remove the previous result file: if GAMS failed on this day, reading a
    % stale result.txt would store the previous day's solution unnoticed.
    if isfile('result.txt'), delete('result.txt'); end
    tic
    [a,b]=system(ejecuta);
    elapsed_time = toc;
    assert(isfile('result.txt'), 'GAMS no genero result.txt en el dia %d', d);

    % Read and store results
    %-------------------------------------------------------------
    % Open the file
    fileID = fopen('result.txt', 'r');
    
    % Read the first line to get the best objective value
    fgetl(fileID); % Skip the first line ('Best objective value achieved')
    bestObjectiveValue = fscanf(fileID, '%f', 1); % Read the numeric value
    
    % Display the best objective value
    fprintf('Best Objective Value on day %i: %.2f\n', d, bestObjectiveValue);
    
    % Read until 'Solution: variable u(i=row, t=col)'
    while ~feof(fileID)
        line = fgetl(fileID);
        if contains(line, 'Solution: variable u(i=row, t=col)')
            break;
        end
    end
    
    % Read matrix u
    u = fscanf(fileID, '%f', [24, 16]); % assuming 16 columns and 24 rows
    u = u';
    
    % Read until 'Solution: variable delta(i=row, t=col)'
    while ~feof(fileID)
        line = fgetl(fileID);
        if contains(line, 'Solution: variable delta(i=row, t=col)')
            break;
        end
    end
    
    % Read matrix delta
    delta = fscanf(fileID, '%f', [24, 16]);
    delta = delta';
    
    % Read until 'Solution: variable P(i=row, t=col)'
    while ~feof(fileID)
        line = fgetl(fileID);
        if contains(line, 'Solution: variable P(i=row, t=col)')
            break;
        end
    end
    
    % Read matrix P
    P = fscanf(fileID, '%f', [24, 16]);
    P = P';
    
    % Read until 'Solution: variable wc(t=col)'
    while ~feof(fileID)
        line = fgetl(fileID);
        if contains(line, 'Solution: variable wc(t=col)')
            break;
        end
    end
    
    % Read vector wc
    wc = fscanf(fileID, '%f', [24, 1]);
    
    % Close the file
    fclose(fileID);
    
    % Store the matrices and best objective value
    Z(d) = bestObjectiveValue;
    U(d,:,:) = u;
    DELTA(d,:,:) = delta;
    POWER(d,:,:) = P;
    WC(d,:) = wc';
    T(d) = elapsed_time;
    
    %-------------------------------------------------------------

    % Load the log file
    fileID = fopen('MIP_GCI_SIMPL_sinPHSP.log', 'r');
    logData = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    
    % Extract the log data as a cell array of strings
    logLines = logData{1};
    
    % Initialize variables to store the elapsed time and relative gap
    elapsedTime = [];
    relativeGap = [];
    
    % Loop through the log lines to find relevant information
    for i = 1:length(logLines)
        line = logLines{i};
        
        % Check for 'elapsed' keyword
        if contains(lower(line), 'elapsed')
            % Extract the numeric value for elapsed time (assuming it's at the end of the line)
            tokens = regexp(line, '(\d+\.?\d*)', 'match');
            if ~isempty(tokens)
                elapsedTime = tokens{end}; % Assume the last number in the line is the elapsed time
            end
        end
        
        % Check for 'relative gap' keyword
        if contains(lower(line), 'relative gap')
            % Extract the numeric value for relative gap (assuming it's at the end of the line)
            tokens = regexp(line, '(\d+\.?\d*)', 'match');
            if ~isempty(tokens)
                relativeGap = tokens{end}; % Assume the last number in the line is the relative gap
            end
        end
    end
    
    % Store the solver statistics
    COMP_TIME(d) = str2double(elapsedTime);
    OPT_GAP(d) = str2double(relativeGap);

    %-------------------------------------------------------------
end

