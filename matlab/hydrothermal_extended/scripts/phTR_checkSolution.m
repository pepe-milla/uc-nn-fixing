function [ok, nViol, why] = phTR_checkSolution(res, NDday)
%PHTR_CHECKSOLUTION  Decides whether a day's solution can be trusted.
%
%   Three failure modes, all of which have actually occurred in this project
%   and none of which announces itself:
%
%   1. GAMS writes result.txt from the put file even when the model was not
%      solved, filling it with the .l values it happens to hold - all zeros
%      for an infeasible instance. Two days of a 730-day campaign were
%      stored that way as Z = 0 with every unit off, dragging the mean cost
%      down unnoticed. A solved day has a positive cost.
%   2. The power balance must close. It is the cheapest way to catch a
%      result.txt that was parsed into the wrong variables.
%   3. The thermal schedule must satisfy the ramps and minimum times. If it
%      does not, GAMS solved a different model - which is exactly what
%      happened when a campaign ran in the folder of the simplified model
%      for 731 days without a single error message.
%
%   The PSH units are not checked: they carry no minimum times and no ramp
%   limits by design.

Pmin = [6.79 6.79 27.84 27.84 75.5 75.5 6.79 6.79 14.09 14.09 4.58 4.58 4.58 6.79 17.7 17.7]';
RU   = [24.26 24.26 29.70 29.70 103.05 113.50 24.26 24.26 20.51 20.51 8.51 8.51 8.51 13.23 22.22 22.22]';
TU   = [1 1 7 7 4 4 1 1 1 1 1 1 1 1 6 6];
TD   = [1 1 7 7 3 3 1 1 1 1 1 1 1 1 5 5];
P0   = [32.3 32.3 37 0 0 0 32.3 32.3 20.5 20.5 8.51 8.51 8.51 17.64 37 0]';
u0   = [1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0]';

ok = false; nViol = NaN; why = '';

if ~(res.z > 0)
    why = sprintf('Z = %.2f no es positivo', res.z);
    return
end

gen = sum(res.P,1) + sum(res.Pg,1) - sum(res.Pp,1);
err = max(abs(gen - (NDday(:)' + res.wc(:)')));
if err > 1
    why = sprintf('balance de potencia descuadrado en %.2f MW', err);
    return
end

u = round(res.u);  delta = round(res.delta);  P = res.P;
uprev = [u0, u(:,1:end-1)];  Pprev = [P0, P(:,1:end-1)];
w  = delta - (u - uprev);

vU = (P - Pprev) - (RU.*uprev + Pmin.*delta);      % RampUP
vD = (Pprev - P) - (RU.*u     + Pmin.*w);          % RampDOWN (RD = RU, SD = Pmin)
nViol = nnz(vU > 0.05) + nnz(vD > 0.05);
nViol = nViol + nnz(u == 1 & P < Pmin - 1e-4) + nnz(u == 0 & P > 1e-4);
for i = 1:16
    for t = 1:24
        if t >= TU(i) && sum(delta(i, max(1,t-TU(i)+1):t)) > u(i,t) + 1e-6
            nViol = nViol + 1;
        end
        if t >= TD(i) && sum(w(i, max(1,t-TD(i)+1):t)) > 1 - u(i,t) + 1e-6
            nViol = nViol + 1;
        end
    end
end

if nViol > 0
    why = sprintf('%d violaciones de rampa o tiempo minimo: GAMS no resolvio este modelo', nViol);
    return
end

ok = true;
end
