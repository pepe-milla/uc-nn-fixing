* Storage of the LP-relaxation solution (Strategy 3 input), thermal case.
* Note the 6-decimal format: the MIP results file writes u with :8:0, which
* would round the relaxed commitments to integers and destroy the very
* information Strategy 3 feeds to the network.

Parameters
rr_z
rr_u(i,t)
;

rr_z = z.l;
rr_u(i,t) = u.l(i,t);

file resultr    /result_rmip.txt/;
put resultr
         put "Relaxed objective value" /;
         put rr_z /;
         put /;

         put "Solution: variable u(i=row, t=col)" /;
         loop(i, put loop(t, put rr_u(i,t):10:6) put /);
