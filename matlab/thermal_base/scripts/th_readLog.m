function [gap, csec] = th_readLog(logFile)
% Extracts the last elapsed time and relative gap reported by CPLEX in the
% GAMS log of the thermal model.
gap = NaN; csec = NaN;
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
end
