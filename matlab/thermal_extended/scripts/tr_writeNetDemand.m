function tr_writeNetDemand(dataGmsFolder, DayDemand, DayWind)
%TR_WRITENETDEMAND  Writes net_demand.gms for one day of the thermal model
%   with intertemporal constraints.
%
%   wcpen = 0, unlike the 1e6 of the simplified thermal case. Curtailment is
%   forced on most days of this formulation - with minimum up/down times the
%   committed units cannot always shut down on windy, low-demand nights - so
%   a penalty of 1e6 EUR/MWh would dominate the objective and make every cost
%   comparison meaningless. At zero the curtailment still carries its real
%   cost, because the thermal fleet has to cover ND + wc, so the model
%   minimises it anyway; what disappears is the artificial spike.
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
assert(fid ~= -1, 'No se pudo escribir net_demand.gms en %s', dataGmsFolder);
fprintf(fid, '* Data reading from Matlab (thermal_extended)\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /0/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end
