function m = sm_parseCplexLog(logFile)
%SM_PARSECPLEXLOG  Solver metrics from a GAMS 23.5 / CPLEX 12.2 log.
%
%   m = SM_PARSECPLEXLOG(logFile) returns a struct with:
%     binsPresolve  binaries left after the first MIP presolve
%     rootLP        objective of the root LP relaxation (first node-log line)
%     rootGap       gap [%] at the end of root processing (last node-0 line
%                   that reports one; NaN if no incumbent existed by then)
%     rootTime      real time of root processing [s]
%     bcTime        real time of branch-and-cut [s] (0 if solved at the root)
%     nodes         nodes explored            (from "MIP Solution: ...")
%     iterations    simplex iterations        (idem)
%     cuts          total cuts applied, and cutsByType (struct)
%     incumbents    number of improved solutions found (node-log lines with *)
%     firstIncNode  node at which the first incumbent was found
%     maxLeft       largest number of open nodes seen in the node log
%     bestBound     "Best possible" at the end
%     relGap        "Relative gap" at the end
%     status        MIP status text (e.g. "integer optimal, tolerance")
%
%   Everything is read from the log only: the models are not modified, and
%   the same function works for the four formulations. Fields that do not
%   appear in a given log (e.g. no branching) are NaN or 0 as documented.

txt = fileread(logFile);
L = regexp(txt, '\r?\n', 'split');

m = struct('binsPresolve',NaN,'rootLP',NaN,'rootGap',NaN,'rootTime',NaN,'bcTime',0, ...
    'nodes',NaN,'iterations',NaN,'cuts',0,'cutsByType',struct(),'incumbents',0, ...
    'firstIncNode',NaN,'maxLeft',0,'bestBound',NaN,'relGap',NaN,'status','');

tok = regexp(txt, 'Reduced MIP has \d+ rows, \d+ columns.*?Reduced MIP has (\d+) binaries', 'tokens', 'once');
if isempty(tok), tok = regexp(txt, 'Reduced MIP has (\d+) binaries', 'tokens', 'once'); end
if ~isempty(tok), m.binsPresolve = str2double(tok{1}); end

tok = regexp(txt, 'Root node processing \(before b&c\):\s*Real time\s*=\s*([\d.]+)', 'tokens', 'once');
if ~isempty(tok), m.rootTime = str2double(tok{1}); end
tok = regexp(txt, 'b&c, \d+ threads?:\s*Real time\s*=\s*([\d.]+)', 'tokens', 'once');
if ~isempty(tok), m.bcTime = str2double(tok{1}); end

tok = regexp(txt, 'MIP Solution:\s*[-\d.eE+]+\s*\((\d+) iterations, (\d+) nodes\)', 'tokens', 'once');
if ~isempty(tok)
    m.iterations = str2double(tok{1});
    m.nodes      = str2double(tok{2});
end
tok = regexp(txt, 'Best possible:\s*([-\d.eE+]+)', 'tokens', 'once');
if ~isempty(tok), m.bestBound = str2double(tok{1}); end
tok = regexp(txt, 'Relative gap:\s*([-\d.eE+]+)', 'tokens', 'once');
if ~isempty(tok), m.relGap = str2double(tok{1}); end
tok = regexp(txt, 'MIP status\(\d+\):\s*([^\r\n]+)', 'tokens', 'once');
if ~isempty(tok), m.status = strtrim(tok{1}); end

% Cuts applied, by type
C = regexp(txt, '([A-Za-z][A-Za-z -]*?) cuts applied:\s*(\d+)', 'tokens');
for k = 1:numel(C)
    name = matlab.lang.makeValidName(strtrim(C{k}{1}));
    n = str2double(C{k}{2});
    if isfield(m.cutsByType, name), m.cutsByType.(name) = m.cutsByType.(name) + n;
    else, m.cutsByType.(name) = n; end
    m.cuts = m.cuts + n;
end

% Node log. A node line is:  [*] Node[+] Left  <objective|integral|infeasible|cutoff> ... [gap%]
first = true;
for k = 1:numel(L)
    s = L{k};
    t = regexp(s, '^\s*(\*?)\s*(\d+)\+?\s+(\d+)\s+(\S+)', 'tokens', 'once');
    if isempty(t), continue; end
    node = str2double(t{2});  left = str2double(t{3});
    if isnan(node) || isnan(left), continue; end
    if first && node == 0
        v = str2double(t{4});
        if ~isnan(v), m.rootLP = v; end
        first = false;
    end
    m.maxLeft = max(m.maxLeft, left);
    if ~isempty(t{1})
        m.incumbents = m.incumbents + 1;
        if isnan(m.firstIncNode), m.firstIncNode = node; end
    end
    if node == 0
        g = regexp(s, '([\d.]+)%\s*$', 'tokens', 'once');
        if ~isempty(g), m.rootGap = str2double(g{1}); end
    end
end
end
