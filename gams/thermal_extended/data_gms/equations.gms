* Data regarding the equations

Equations

*Equations concerning the TGUs
MaxPowerTGU(i,q,t)      maximum output power of TGU i corresponding to the qth segment of the production cost curve
LogicTGU1(i,t)          logic of the TGU binary variables (1)
LogicTGU2(i,t)          logic of the TGU binary variables (2)
PowerVarTGU1(i,t)       calculates the power variation of TGU i between period t-1 and t
PowerVarTGU2(i,t)       calculates the power variation of TGU i between period t-1 and t
MinTimeOnTGU(i,t)       minimum time during which unit i must remain on-line after starting-up
MinTimeOffTGU(i,t)      minimum time during which unit i must remain off-line after shutting-down
InitialState(i,t)       minimum time during which the unit i must remain on or off at the beggining of the scheduling horizon
RampUP(i,t)             ramp constraints up
RampDOWN(i,t)           ramp constrants down

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
LogicTGU1(i,t).. delta(i,t) - w(i,t) =e= u(i,t) - u(i,t-1)$(ord(t) ne 1) - u0(i)$(ord(t) eq 1);
LogicTGU2(i,t).. delta(i,t) + w(i,t) =l= 1;

* Max output power of each segment of the production cost curve of each TGU
MaxPowerTGU(i,q,t).. p(i,q,t) =l= u(i,t)*pinc(i,q);

* Minimum time on/off after a TGU unit starts-up/shuts-down
MinTimeOnTGU(i,t)$(ord(t) ge TU(i)).. sum(tt$(ord(tt) gt ord(t)-TU(i) and ord(tt) le ord(t)),delta(i,tt)) =l= u(i,t);
MinTimeOffTGU(i,t)$(ord(t) ge TD(i)).. sum(tt$(ord(tt) gt ord(t)-TD(i) and ord(tt) le ord(t)),w(i,tt)) =l= 1 - u(i,t);
InitialState(i,t)$(TUR(i)+TDR(i) ge 1 and ord(t) le TUR(i)+TDR(i)).. u(i,t) =e= u0(i);

* Calculation of the power variation of the TGUs. Note: Here "=g=" is used instead of "=e=", but considering that pv is in the cost function, the model will try to adjust it to its lower bound.
PowerVarTGU1(i,t).. pv(i,t) =g= (sum(q,p(i,q,t))-sum(q,p(i,q,t-1)))$(ord(t) ne 1) + (sum(q,p(i,q,t))-(P0(i)-Pmin(i)$(P0(i) ne 0)))$(ord(t) eq 1);
PowerVarTGU2(i,t).. pv(i,t) =g= -(sum(q,p(i,q,t))-sum(q,p(i,q,t-1)))$(ord(t) ne 1) - (sum(q,p(i,q,t))-(P0(i)-Pmin(i)$(P0(i) ne 0)))$(ord(t) eq 1);

* Ramping constraints
RampUP(i,t)..    u(i,t)*Pmin(i) + sum(q,p(i,q,t))
               - (u(i,t-1)*Pmin(i) + sum(q,p(i,q,t-1)))$(ord(t) ne 1)
               - P0(i)$(ord(t) eq 1)
               =l= RU(i)*(u(i,t-1)$(ord(t) ne 1) + u0(i)$(ord(t) eq 1))
                 + SU(i)*delta(i,t);

RampDOWN(i,t)..  (u(i,t-1)*Pmin(i) + sum(q,p(i,q,t-1)))$(ord(t) ne 1)
               + P0(i)$(ord(t) eq 1)
               - (u(i,t)*Pmin(i) + sum(q,p(i,q,t)))
               =l= RD(i)*u(i,t) + SD(i)*w(i,t);

* Spinning reserve regarding the expected inter-hourly increase in power demand and the most likely wind power loss
SpReserve01(t).. sum(i,u(i,t)*sum(q,pinc(i,q))-sum(q,p(i,q,t))) =g= max(sr_d(t),sr_wp(t));

* Spinning reserve regarding the loss of the largest TGU
SpReserve02(i,t).. sum(ii$(ord(ii) ne ord(i)),u(ii,t)*sum(q,pinc(ii,q))-sum(q,p(ii,q,t)))=g= u(i,t)*Pmin(i)+sum(q,p(i,q,t));

