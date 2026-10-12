
$Title SIMPLIFIED MIP Model used in Great Canary (DAUC) without the PHES station of Chira-Soria

$include /data_gms/set.gms
$include /data_gms/net_demand.gms
$include /data_gms/spinning_reserve.gms
$include /data_gms/data_TGUs.gms
$include /data_gms/initial_conditions.gms
$include /data_gms/variables.gms
$include /data_gms/fixed_variables.gms
$include /data_gms/equations.gms

Model UC uses all constraints /all/;
UC.optfile=9;

Solve UC using MIP minimizing z;

$include /data_gms/results.gms

