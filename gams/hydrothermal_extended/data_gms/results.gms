* Code for the analysis and storage of results

Parameters
r_z
r_u(i,t)
r_ug(c,t)
r_up(c,t)
r_delta(i,t)
r_deltag(c,t)
r_deltap(c,t)
r_P(i,t)
r_Pg(c,t)
r_Pp(c,t)
r_v(t)
r_wc(t)
;

r_z = z.l;
r_u(i,t) = u.l(i,t);
r_ug(c,t) = ug.l(c,t);
r_up(c,t) = up.l(c,t);
r_delta(i,t) = delta.l(i,t);
r_deltag(c,t) = deltag.l(c,t);
r_deltap(c,t) = deltap.l(c,t);
r_P(i,t)  = u.l(i,t)*Pmin(i) + sum(q,p.l(i,q,t));
r_Pg(c,t) = ug.l(c,t)*Pgmin(c) + pg.l(c,t);
r_Pp(c,t) = up.l(c,t)*Ppmin(c) + pp.l(c,t);
r_v(t) = v.l(t);
r_wc(t) = wc.l(t);

file result    /result.txt/;
put result
         put "Best objetive value achieved" /;
         put r_z /;
         put /;

         put "Solution: variable u(i=row, t=col)" /;
         loop(i, put loop(t, put r_u(i,t):8:4) put /);
         put /;

         put "Solution: variable ug(c=row, t=col)" /;
         loop(c, put loop(t, put r_ug(c,t):8:4) put /);
         put /;

         put "Solution: variable up(c=row, t=col)" /;
         loop(c, put loop(t, put r_up(c,t):8:4) put /);
         put /;

         put "Solution: variable delta(i=row, t=col)" /;
         loop(i, put loop(t, put r_delta(i,t):8:4) put /);
         put /;

         put "Solution: variable deltag(c=row, t=col)" /;
         loop(c, put loop(t, put r_deltag(c,t):8:4) put /);
         put /;

         put "Solution: variable deltap(c=row, t=col)" /;
         loop(c, put loop(t, put r_deltap(c,t):8:4) put /);
         put /;

         put "Solution: variable P(i=row, t=col)" /;
         loop(i, put loop(t, put r_P(i,t):8:2) put /);
         put /;

         put "Solution: variable Pg(c=row, t=col)" /;
         loop(c, put loop(t, put r_Pg(c,t):8:2) put /);
         put /;

         put "Solution: variable Pp(c=row, t=col)" /;
         loop(c, put loop(t, put r_Pp(c,t):8:2) put /);
         put /;

         put "Solution: variable v(t=col)" /;
         loop(t, put r_v(t):8:0);
         put /;

         put "Solution: variable wc(t=col)" /;
         loop(t, put r_wc(t):8:2);
