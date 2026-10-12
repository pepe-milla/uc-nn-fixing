function th_writeNetDemand(dataGmsFolder, DayDemand, DayWind)
% Writes net_demand.gms for one day, with the same format as
% write_net_demand.m but into an explicit folder (so that each parallel
% worker writes into its own copy of the model).
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
fprintf(fid, '* Data reading from Matlab (thermal case, parallel campaign)\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /1000000/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end
