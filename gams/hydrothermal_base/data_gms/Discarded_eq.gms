Equations

MinTimeOn(i,t)          minimum time during which TGU i must remain on
MinTimeOff(i,t)         minimum time during which TGU i must remain off
;
* Minimum time on and off of each TGU
MinTimeOn(i,t)$(ord(t) ge TU(i)).. sum(tt$((ord(tt)ge(ord(t)-TU(i)+1))and(ord(tt)le(ord(t)))),delta(i,tt)) =l= u(i,t);

Parameter

TU(i) minimum time during which the unit must remain on (h)
/
1 1
2 1
3 7
4 7
5 4
6 4
7 1
8 1
9 1
10 1
11 1
12 1
13 1
14 1
15 6
16 6
/
TD(i) minimum time during which the unit must remain off (h)
/
1 1
2 1
3 7
4 7
5 3
6 3
7 1
8 1
9 1
10 1
11 1
12 1
13 1
14 1
15 5
16 5
/
TU0(i) Number of hours that the unit has been up before the scheduling horizon (h)
/
1 10
2 10
3 10
4 0
5 0
6 0
7 10
8 10
9 10
10 10
11 10
12 10
13 10
14 10
15 10
16 0
17 0
18 0
19 0
20 0
21 0
22 0
23 0
24 0
25 0
26 0
27 0
28 0
/
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
17 0
18 0
19 0
20 0
21 0
22 0
23 0
24 0
25 0
26 0
27 0
28 0
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
17 0
18 0
19 0
20 0
21 0
22 0
23 0
24 0
25 0
26 0
27 0
28 0
/
TD0(i) Number of hours that the unit has been down before the scheduling horizon (h)
/
1 0
2 0
3 0
4 10
5 10
6 240
7 0
8 0
9 0
10 0
11 0
12 0
13 0
14 0
15 0
16 10
17 240
18 240
19 240
20 240
21 240
22 240
23 240
24 240
25 240
26 240
27 240
28 240
/
;