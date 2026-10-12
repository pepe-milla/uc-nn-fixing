* Data regarding the equations

Equations

*Equations concerning the TGUs
MaxPowerTGU(i,q,t)      maximum output power of TGU i corresponding to the qth segment of the production cost curve
LogicTGU(i,t)           logic of the TGU binary variables
PowerVarTGU1(i,t)       calculates the power variation of TGU i between period t-1 and t
PowerVarTGU2(i,t)       calculates the power variation of TGU i between period t-1 and t

*Equations concerning all units
TotalCost               aggregate power supply cost (including penalty for wind curtailments)
PowerBal(t)             total power balance
SpReserve01(t)          spinning reserve regarding demand & windpower
SpReserve02(i,t)        spinning reserve regarding a generator failure
;

* Objective function
*-------------------------------------------------------------------

TotalCost..
z =e= sum(t, sum(i, u(i,t)*Cf(i)+sum(q,p(i,q,t)*Cv(i,q)) + delta(i,t)*Mf(i)+(u(i,t)*Pmin(i)+sum(q,p(i,q,t)))*Mv(i) + delta(i,t)*Csu(i)+pv(i,t)*Pvc(i)) + wcpen*wc(t));


* Model restrictions (in the same order as in Fernández-Muñoz(2020b))
*-------------------------------------------------------------------

* Power balance
PowerBal(t).. D(t)+wc(t) =e= sum(i,u(i,t)*Pmin(i)+sum(q,p(i,q,t)));

* Logic constraint on the startup binary variable for the TGUs
LogicTGU(i,t).. delta(i,t) =g= u(i,t) - u(i,t-1)$(ord(t) ne 1) - u0(i)$(ord(t) eq 1);

* Max output power of each segment of the production cost curve of each TGU
MaxPowerTGU(i,q,t).. p(i,q,t) =l= u(i,t)*pinc(i,q);

* Calculation of the power variation of the TGUs. Note: Here "=g=" is used instead of "=e=", but considering that pv is in the cost function, the model will try to adjust it to its lower bound.
PowerVarTGU1(i,t).. pv(i,t) =g= (sum(q,p(i,q,t))-sum(q,p(i,q,t-1)))$(ord(t) ne 1) + (sum(q,p(i,q,t))-(P0(i)-Pmin(i)$(P0(i) ne 0)))$(ord(t) eq 1);
PowerVarTGU2(i,t).. pv(i,t) =g= -(sum(q,p(i,q,t))-sum(q,p(i,q,t-1)))$(ord(t) ne 1) - (sum(q,p(i,q,t))-(P0(i)-Pmin(i)$(P0(i) ne 0)))$(ord(t) eq 1);

* Spinning reserve regarding the expected inter-hourly increase in power demand and the most likely wind power loss
SpReserve01(t).. sum(i,u(i,t)*sum(q,pinc(i,q))-sum(q,p(i,q,t))) =g= max(sr_d(t),sr_wp(t));

* Spinning reserve regarding the loss of the largest TGU
SpReserve02(i,t).. sum(ii$(ord(ii) ne ord(i)),u(ii,t)*sum(q,pinc(ii,q))-sum(q,p(ii,q,t)))=g= u(i,t)*Pmin(i)+sum(q,p(i,q,t));

