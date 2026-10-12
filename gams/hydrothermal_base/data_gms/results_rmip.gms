* Storage of the LP-relaxation solution (Strategy 3 input)

Parameters
rr_z
rr_u(i,t)
rr_ug(c,t)
rr_up(c,t)
;

rr_z = z.l;
rr_u(i,t)  = u.l(i,t);
rr_ug(c,t) = ug.l(c,t);
rr_up(c,t) = up.l(c,t);

file resultr    /result_rmip.txt/;
put resultr
         put "Relaxed objective value" /;
         put rr_z /;
         put /;

         put "Solution: variable u(i=row, t=col)" /;
         loop(i, put loop(t, put rr_u(i,t):10:6) put /);
         put /;

         put "Solution: variable ug(c=row, t=col)" /;
         loop(c, put loop(t, put rr_ug(c,t):10:6) put /);
         put /;

         put "Solution: variable up(c=row, t=col)" /;
         loop(c, put loop(t, put rr_up(c,t):10:6) put /);
