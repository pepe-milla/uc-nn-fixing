
% Initial and final day of the training period
init_time = datetime('2021-12-19 20:00','InputFormat','yyyy-MM-dd HH:mm');
final_time = datetime('2023-12-19 20:00','InputFormat','yyyy-MM-dd HH:mm');

% Baseline guard: the model includes data_gms/fixed_variables.gms
% unconditionally, so any content left there by a previous NN campaign
% would silently constrain every instance solved here. Emptied once,
% before the loop.
modelDir = fullfile(uc_config().gamsRoot, 'thermal_extended');
fid = fopen(fullfile(modelDir,'data_gms','fixed_variables.gms'), 'w');
fprintf(fid, '* No fixed variables (label generation)\n');
fclose(fid);
[DAYDEM, DAYWIND] = uc_load_profiles('2021-2023');
numDays = days(duration(final_time-init_time));
DayDemand = zeros(1,24);
DayWind = zeros(1,24);

% Result variables for storage. Preallocated with NaN, not zeros: a day that
% GAMS could not solve must stand out instead of passing for a day whose
% cost happened to be zero. That is how days 660 and 661 of the 2021-2023
% campaign went unnoticed - both infeasible, both stored as Z = 0 with every
% unit off, and both dragging the mean cost down.
Z = nan(numDays,1);
U = nan(numDays,16,24);
DELTA = nan(numDays,16,24);
POWER = nan(numDays,16,24);
WC = nan(numDays,24);
COMP_TIME = nan(numDays,1);
OPT_GAP = nan(numDays,1);
T = nan(numDays,1);
ND = zeros(numDays,24);
INFEAS = false(numDays,1);


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
    T(d) = elapsed_time;

    % Solver statistics and failure detection. GAMS writes result.txt from
    % the put file even when the model was not solved, filling it with the
    % .l values it happens to hold - all zeros for an infeasible instance -
    % so the existence of the file proves nothing.
    %
    % Do NOT match the bare word "infeasible": the CPLEX node log prints it
    % in the status column of any node whose LP relaxation is infeasible,
    % which is routine and would discard perfectly good days. The phrasings
    % below appear only in the solver's verdict. "*** Status: Normal
    % completion" is no use either - it shows up even for an infeasible
    % model, since it reports that GAMS itself finished.
    logLines = strings(0,1);
    if isfile('MIP_GCI_SIMPL_sinPHSP.log')
        logLines = readlines('MIP_GCI_SIMPL_sinPHSP.log');
    end
    INFEAS(d) = any(contains(lower(logLines), 'integer infeasible')) || ...
                any(contains(lower(logLines), 'problem is infeasible'));

    % Elapsed time and relative gap. GAMS writes the job time as H:MM:SS.sss,
    % so the fields must be combined: reading only the last number of the
    % line, as an earlier version of this parser did, silently dropped the
    % minutes and reported a 2 min 45 s solve as 45 s. Everything under a
    % minute came out right, which is what made it hard to notice. The
    % COMP_TIME stored by the campaigns run before this fix is capped just
    % under 60 s and must not be used; T, measured with tic/toc, is fine.
    OPT_GAP(d) = NaN; COMP_TIME(d) = NaN;
    tok = regexp(strjoin(cellstr(logLines), newline), 'elapsed\s+([\d:]+\.\d+)', 'tokens');
    if ~isempty(tok)
        parts = sscanf(tok{end}{1}, '%d:%d:%f');
        COMP_TIME(d) = parts(1)*3600 + parts(2)*60 + parts(3);
    end
    tok = regexp(strjoin(cellstr(logLines), newline), 'Relative gap:\s+([\d\.]+)', 'tokens');
    if ~isempty(tok)
        OPT_GAP(d) = str2double(tok{end}{1});
    end

    if INFEAS(d) || ~isfile('result.txt')
        warning('Dia %d: GAMS no devolvio una solucion utilizable (INFEAS = %d). Se deja en NaN.', ...
            d, INFEAS(d));
        continue
    end

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
    
    % Last line of defence: a solved day always has a positive cost and
    % covers its demand. This catches the case where GAMS reports normal
    % completion but the put file wrote an empty solution.
    if ~(bestObjectiveValue > 0) || max(abs(sum(P,1) - (ND(d,:) + wc'))) > 1
        warning('Dia %d: solucion incoherente (Z = %.2f, balance roto). Se deja en NaN.', ...
            d, bestObjectiveValue);
        INFEAS(d) = true;
        continue
    end

    % Store the matrices and best objective value
    Z(d) = bestObjectiveValue;
    U(d,:,:) = u;
    DELTA(d,:,:) = delta;
    POWER(d,:,:) = P;
    WC(d,:) = wc';

    %-------------------------------------------------------------
end

fprintf('\n%d dias, %d sin solucion utilizable.\n', numDays, nnz(INFEAS));
if any(INFEAS)
    fprintf('Dias afectados: %s\n', mat2str(find(INFEAS)'));
    fprintf('Quedan como NaN: excluyelos de las medias o resuelvelos aparte.\n');
end
fprintf('Coste medio %.0f EUR, tiempo medio %.2f s (sobre los dias validos).\n', ...
    mean(Z,'omitnan'), mean(COMP_TIME,'omitnan'));

