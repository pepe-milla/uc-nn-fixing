
$Title LP RELAXATION of the SIMPLIFIED MIP Model used in Great Canary (DAUC) without the PHES station of Chira-Soria

* Same model as MIP_GCI_SIMPL_sinPHSP.gms with two differences:
*   - fixed_variables.gms is NOT included: the relaxation is solved on the
*     unrestricted problem, since its solution is an INPUT to the network
*   - it is solved as RMIP, so all binary variables are relaxed to [0,1]
* Results go to result_rmip.txt, leaving result.txt untouched, so that both
* solves can coexist in the same model (or worker) folder.

$include /data_gms/set.gms
$include /data_gms/net_demand.gms
$include /data_gms/spinning_reserve.gms
$include /data_gms/data_TGUs.gms
$include /data_gms/initial_conditions.gms
$include /data_gms/variables.gms
$include /data_gms/equations.gms

Model UC uses all constraints /all/;
UC.optfile=9;

Solve UC using RMIP minimizing z;

$include /data_gms/results_rmip.gms
