function [gap, csec, infeas] = phTR_readLog(logFile)
%PHTR_READLOG  Relative gap, elapsed time and infeasibility from a GAMS log.
%
%   Elapsed time: GAMS writes the job time as H:MM:SS.sss, so the fields
%   have to be combined. Taking only the last number of the line silently
%   drops the minutes and reports a 2 min 45 s solve as 45 s - everything
%   under a minute comes out right, which is what makes it hard to notice.
%
%   Infeasibility: only the solver's verdict counts. The bare word
%   "infeasible" also appears in the status column of any branch-and-bound
%   node whose relaxation is infeasible, which is routine and would discard
%   perfectly good days. Nor is "*** Status: Normal completion" of any use:
%   it reports that GAMS finished, and shows up for an infeasible model too.
gap = NaN; csec = NaN; infeas = false;
if ~isfile(logFile), return; end
txt = fileread(logFile);

tok = regexp(txt, 'elapsed\s+([\d:]+\.\d+)', 'tokens');
if ~isempty(tok)
    parts = sscanf(tok{end}{1}, '%d:%d:%f');
    csec = parts(1)*3600 + parts(2)*60 + parts(3);
end

tok = regexp(txt, 'Relative gap:\s+([\d\.]+)', 'tokens');
if ~isempty(tok)
    gap = str2double(tok{end}{1});
end

low = lower(txt);
infeas = contains(low, 'integer infeasible') || contains(low, 'problem is infeasible');
end
