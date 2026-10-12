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

*Equations concerning the PSHP
MaxPowerTurb(c,t)       maximum output power of PSH unit c in turbine mode
MaxPowerPump(c,t)       maximum input power of  PSH unit c in pump mode
LogicTurb(c,t)          logic of the PSH units in turbine mode
LogicPump(c,t)          logic of the PSH units in pump mode
PowerVarTurb1(c,t)      calculates the power variation of PSH unit c between period t-1 and t in turbine mode
PowerVarTurb2(c,t)      calculates the power variation of PSH unit c between period t-1 and t in turbine mode
PowerVarPump1(c,t)      calculates the power variation of PSH unit c between period t-1 and t in pump mode
PowerVarPump2(c,t)      calculates the power variation of PSH unit c between period t-1 and t in pump mode
WaterBal(t)             balance equation of the upper reservoir
FinalVol(t)             final volume of the uper reservoir
LogicPSHP(c,t)          logic of the PSH binary units to avoid hydraulic short-circuit
*SymTurb(c,t)            symmetry breaking among identical PSH units in turbine mode (stacked order)
*SymPump(c,t)            symmetry breaking among identical PSH units in pump mode (stacked order)

*Equations concerning all units
TotalCost               aggregate power supply cost (including penalty for wind curtailments)
PowerBal(t)             total power balance
SpReserve01(t)          spinning reserve regarding demand & windpower
SpReserve02(i,t)        spinning reserve regarding a generator failure
SpReserve03(c,t)        spinning reserve regarding a generator failure
;

* Objective function
*-------------------------------------------------------------------

TotalCost..
z =e= sum(t, sum(i, u(i,t)*Cf(i)+sum(q,p(i,q,t)*Cv(i,q)) + delta(i,t)*Mf(i)+(u(i,t)*Pmin(i)+sum(q,p(i,q,t)))*Mv(i) + delta(i,t)*Csu(i)+pv(i,t)*Pvc(i))
 + sum(c, deltag(c,t)*Csug(c)+deltap(c,t)*Csup(c)+pvp(c,t)*Pvcp(c)) + (v0$(ord(t) eq 1)+v(t-1)$(ord(t) ne 1)-v(t))*wv + wcpen*wc(t));


* Model restrictions (in the same order as in Fern�ndez-Mu�oz(2020b))
*-------------------------------------------------------------------

* Power balance
PowerBal(t).. D(t)+wc(t) =e= sum(i,u(i,t)*Pmin(i)+sum(q,p(i,q,t))) + sum(c,ug(c,t)*Pgmin(c)+pg(c,t)) - sum(c,up(c,t)*Ppmin(c)+pp(c,t));

* Logic constraint on the startup binary variable for the TGUs. The exact
* equality with the shut-down variable replaces the inequality used in the
* simplified formulation: the minimum-time and ramp constraints below need
* w(i,t) to identify shut-downs, which delta alone cannot do.
LogicTGU1(i,t).. delta(i,t) - w(i,t) =e= u(i,t) - u(i,t-1)$(ord(t) ne 1) - u0(i)$(ord(t) eq 1);
LogicTGU2(i,t).. delta(i,t) + w(i,t) =l= 1;

* Minimum time on/off after a TGU unit starts-up/shuts-down. The PSH units
* are deliberately left out: they have no minimum times.
MinTimeOnTGU(i,t)$(ord(t) ge TU(i)).. sum(tt$(ord(tt) gt ord(t)-TU(i) and ord(tt) le ord(t)),delta(i,tt)) =l= u(i,t);
MinTimeOffTGU(i,t)$(ord(t) ge TD(i)).. sum(tt$(ord(tt) gt ord(t)-TD(i) and ord(tt) le ord(t)),w(i,tt)) =l= 1 - u(i,t);
InitialState(i,t)$(TUR(i)+TDR(i) ge 1 and ord(t) le TUR(i)+TDR(i)).. u(i,t) =e= u0(i);

* Ramping constraints of the TGUs. Again only for the thermal units: the
* pump-turbines reach full output within minutes, so an hourly ramp limit
* would never bind.
RampUP(i,t)..    u(i,t)*Pmin(i) + sum(q,p(i,q,t))
               - (u(i,t-1)*Pmin(i) + sum(q,p(i,q,t-1)))$(ord(t) ne 1)
               - P0(i)$(ord(t) eq 1)
               =l= RU(i)*(u(i,t-1)$(ord(t) ne 1) + u0(i)$(ord(t) eq 1))
                 + SU(i)*delta(i,t);

RampDOWN(i,t)..  (u(i,t-1)*Pmin(i) + sum(q,p(i,q,t-1)))$(ord(t) ne 1)
               + P0(i)$(ord(t) eq 1)
               - (u(i,t)*Pmin(i) + sum(q,p(i,q,t)))
               =l= RD(i)*u(i,t) + SD(i)*w(i,t);

* Logic constraint on the startup binary variables for the PSH units in turbine mode
LogicTurb(c,t).. deltag(c,t) =g= ug(c,t) - ug(c,t-1)$(ord(t) ne 1) - ug0(c)$(ord(t) eq 1);

* Logic constraint on the startup binary variables for the PSH units in pump mode
LogicPump(c,t).. deltap(c,t) =g= up(c,t) - up(c,t-1)$(ord(t) ne 1) - up0(c)$(ord(t) eq 1);

* Max output power of each segment of the production cost curve of each TGU
MaxPowerTGU(i,q,t).. p(i,q,t) =l= u(i,t)*pinc(i,q);

* Max output power of each PSH unit in turbine mode
MaxPowerTurb(c,t).. pg(c,t) =l= ug(c,t)*(Pgmax(c)-Pgmin(c));

* Max output power of each PSH unit in pump mode
MaxPowerPump(c,t).. pp(c,t) =l= up(c,t)*(Ppmax(c)-Ppmin(c));

* Calculation of the power variation of the TGUs. Note: Here "=g=" is used instead of "=e=", but considering that pv is in the cost function, the model will try to adjust it to its lower bound.
PowerVarTGU1(i,t).. pv(i,t) =g= (sum(q,p(i,q,t))-sum(q,p(i,q,t-1)))$(ord(t) ne 1) + (sum(q,p(i,q,t))-(P0(i)-Pmin(i)$(P0(i) ne 0)))$(ord(t) eq 1);
PowerVarTGU2(i,t).. pv(i,t) =g= -(sum(q,p(i,q,t))-sum(q,p(i,q,t-1)))$(ord(t) ne 1) - (sum(q,p(i,q,t))-(P0(i)-Pmin(i)$(P0(i) ne 0)))$(ord(t) eq 1);

* Calculation of the power variation of the PSH units in turbine mode. Note: Here "=g=" is used instead of "=e=", but considering that pv is in the cost function, the model will try to adjust it to its lower bound.
PowerVarTurb1(c,t).. pvp(c,t) =g= (pg(c,t)-pg(c,t-1))$(ord(t) ne 1) + (pg(c,t)-(Pg0(c)-Pgmin(c)$(Pg0(c) ne 0)))$(ord(t) eq 1);
PowerVarTurb2(c,t).. pvp(c,t) =g= -(pg(c,t)-pg(c,t-1))$(ord(t) ne 1) - (pg(c,t)-(Pg0(c)-Pgmin(c)$(Pg0(c) ne 0)))$(ord(t) eq 1);

* Calculation of the power variation of the PSH units in pump mode. Note: Here "=g=" is used instead of "=e=", but considering that pv is in the cost function, the model will try to adjust it to its lower bound.
PowerVarPump1(c,t).. pvp(c,t) =g= (Pp(c,t)-Pp(c,t-1))$(ord(t) ne 1) + (pp(c,t)-(Pp0(c)-Ppmin(c)$(Pp0(c) ne 0)))$(ord(t) eq 1);
PowerVarPump2(c,t).. pvp(c,t) =g= -(Pp(c,t)-Pp(c,t-1))$(ord(t) ne 1) - (pp(c,t)-(Pp0(c)-Ppmin(c)$(Pp0(c) ne 0)))$(ord(t) eq 1);

* Spinning reserve regarding the expected inter-hourly increase in power demand and the most likely wind power loss
SpReserve01(t).. sum(i,u(i,t)*sum(q,pinc(i,q))-sum(q,p(i,q,t))) + sum(c,ug(c,t)*(Pgmax(c)-Pgmin(c))-pg(c,t)) + sum(c,pp(c,t)) =g= max(sr_d(t),sr_wp(t));

* Spinning reserve regarding the loss of the largest TGU
SpReserve02(i,t).. sum(ii$(ord(ii) ne ord(i)),u(ii,t)*sum(q,pinc(ii,q))-sum(q,p(ii,q,t))) + sum(c,ug(c,t)*(Pgmax(c)-Pgmin(c))-pg(c,t)) + sum(c,pp(c,t)) =g= u(i,t)*Pmin(i)+sum(q,p(i,q,t));

* Spinning reserve regarding the loss of the largest PSH turbine unit
SpReserve03(c,t).. sum(i,u(i,t)*sum(q,pinc(i,q))-sum(q,p(i,q,t))) + sum(cc$(ord(cc) ne ord(c)),ug(cc,t)*(Pgmax(cc)-Pgmin(cc))-pg(cc,t)) + sum(cc,pp(cc,t)) =g= ug(c,t)*Pgmin(c)+pg(c,t);

* Water volume balance in the upper reservoir
WaterBal(t).. v(t) =e= v0$(ord(t) eq 1) + v(t-1)$(ord(t) ne 1) + 3600*(sum(c,up(c,t)*qminp(c)*qpn(c)+rp(c)*pp(c,t))-sum(c,ug(c,t)*qmint(c)*qtn(c)+rt(c)*pg(c,t)));

* Ensure that the upper reservoir volume at the end of the period is greater than vf
FinalVol(t)$(ord(t) eq card(t)).. v(t) =g= vf;

* Impossibility of simultaneous modes of operation for the PSH units
LogicPSHP(c,t).. ug(c,t) + up(c,t) =l= 1;

* Symmetry breaking
*SymTurb(c,t)$(ord(c) lt card(c)).. ug(c,t) =g= ug(c+1,t);
*SymPump(c,t)$(ord(c) lt card(c)).. up(c,t) =l= up(c+1,t);

