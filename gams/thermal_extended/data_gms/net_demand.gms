* Sample day. Overwritten by the MATLAB scripts (see matlab/thermal_extended/scripts)
* before each GAMS call.
Scalar
wcpen wind curtailment penalty (EUR:MWh) /0/;
Parameter
Demand(t) Hourly REAL power demands
/
1 487.575
2 464.6167
3 404.3583
4 350.4917
5 316.45
6 294.5083
7 287.975
8 282.1417
9 282.9167
10 297.4333
11 342.3333
12 389.3417
13 408.05
14 431.3917
15 437.5333
16 437.1917
17 439.975
18 445.7417
19 448.9
20 436.1667
21 425.3333
22 433.325
23 464.0583
24 485.2667
/
Wind(t) Scenarios of wind power production
/
1  13.3417
2  16.2667
3  26.4
4  36.025
5  38.7667
6  42.5083
7  36.7
8  34.1667
9  25.275
10  20.35
11  14.5
12  9.0667
13  6.8
14  7.2917
15  13.7
16  23.5083
17  32.6333
18  46.0583
19  72.725
20  63.5667
21  63.4667
22  68.2333
23  40.425
24  54.625
/
;
Parameter
D(t) Net power demand (MW);
D(t)=Demand(t)-Wind(t)
;
