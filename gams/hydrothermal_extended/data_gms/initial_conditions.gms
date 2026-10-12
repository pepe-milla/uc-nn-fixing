*Initial conditions of the TGUs
*calculadas para el d�a 264 de 2020 con horizonte 2025
Scalar
v0 volume of the upper reservoir at initial stage /2000000/
vf volume of the upper reservoir at final stage /2000000/
wv marginal water value  (�:m3) /0/
*0.085
;
Parameter
P0(i) Power dispatched  at the beginning of time horizon by the TGUs [MW]
/
1 32.3
2 32.3
3 37
4 0
5 0
6 0
7 32.3
8 32.3
9 20.5
10 20.5
11 8.51
12 8.51
13 8.51
14 17.64
15 37
16 0
/
Pg0(c) Power dispatched  at the beginning of time horizon by the PSH units in turbine mode [MW]
/
1 0
2 0
3 0
4 0
5 0
6 0
/
Pp0(c) Power consumed  at the beginning of time horizon by the PSH units in pump mode [MW]
/
1 0
2 0
3 0
4 0
5 0
6 0
/
u0(i) Initial unit status of TGUs (1 = 0n - 0 = off)
/
1 1
2 1
3 1
4 0
5 0
6 0
7 1
8 1
9 1
10 1
11 1
12 1
13 1
14 1
15 1
16 0
/
ug0(c) Initial unit status of PSH units in turbine mode (1 = 0n - 0 = off)
/
1 0
2 0
3 0
4 0
5 0
6 0
/
up0(c) Initial unit status of PSH units in pump mode (1 = 0n - 0 = off)
/
1 0
2 0
3 0
4 0
5 0
6 0
/
* Residual obligations inherited from the previous day. Kept at zero, as in
* "Modelo GC simpl_T_Ramp_sinPSHP": the initial state is synthetic and the
* same for every day, so carrying an obligation over would be arbitrary.
TUR(i) number of initial hours during which the unit must remain up at the beginning of the scheduling horizon (h)
/
1 0
2 0
3 0
4 0
5 0
6 0
7 0
8 0
9 0
10 0
11 0
12 0
13 0
14 0
15 0
16 0
/
TDR(i) number of initial hours during which the unit must remain down at the beginning of the scheduling horizon (h)
/
1 0
2 0
3 0
4 0
5 0
6 0
7 0
8 0
9 0
10 0
11 0
12 0
13 0
14 0
15 0
16 0
/
;
