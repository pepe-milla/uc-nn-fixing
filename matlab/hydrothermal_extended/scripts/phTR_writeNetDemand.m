function phTR_writeNetDemand(dataGmsFolder, DayDemand, DayWind)
%PHTR_WRITENETDEMAND  Writes net_demand.gms for one day of the hydrothermal
%   model with intertemporal constraints.
%
%   wcpen = 0. The penalty never activates in the hydrothermal model - the
%   pumps absorb any surplus, and there was zero curtailment in the 3650
%   day-solutions of the stored PHES campaigns - so the value does not
%   change any cost observed so far and the results stay comparable with
%   those runs. It is zero as insurance: pumping is capped by the volume of
%   the upper reservoir, and with minimum times forcing the slow units to
%   keep running overnight the reservoir may fill sooner. The day
%   curtailment does become unavoidable, a penalty of 1e6 EUR/MWh would turn
%   that day's cost into millions and wreck every average of the campaign.
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
assert(fid ~= -1, 'No se pudo escribir net_demand.gms en %s', dataGmsFolder);
fprintf(fid, '* Data reading from Matlab (hydrothermal_extended)\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /0/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end
