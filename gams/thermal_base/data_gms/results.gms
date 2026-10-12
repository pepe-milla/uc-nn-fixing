* Code for the analysis and storage of results

Parameters
r_z
r_u(i,t)
r_delta(i,t)
r_P(i,t)
r_wc(t)
;

r_z = z.l;
r_u(i,t) = u.l(i,t);
r_delta(i,t) = delta.l(i,t);
r_P(i,t)  = u.l(i,t)*Pmin(i) + sum(q,p.l(i,q,t));
r_wc(t) = wc.l(t);

file result    /result.txt/;
put result
         put "Best objetive value achieved" /;
         put r_z /;
         put /;

         put "Solution: variable u(i=row, t=col)" /;
         loop(i, put loop(t, put r_u(i,t):8:0) put /);
         put /;

         put "Solution: variable delta(i=row, t=col)" /;
         loop(i, put loop(t, put r_delta(i,t):8:0) put /);
         put /;

         put "Solution: variable P(i=row, t=col)" /;
         loop(i, put loop(t, put r_P(i,t):8:2) put /);
         put /;

         put "Solution: variable wc(t=col)" /;
         loop(t, put r_wc(t):8:2);
