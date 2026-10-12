* Sample day. Overwritten by the MATLAB scripts (see matlab/hydrothermal_base/scripts)
* before each GAMS call.
Scalar
wcpen wind curtailment penalty (EUR:MWh) /1000000/;
Parameter
Demand(t) Hourly REAL power demands
/
1 455.75
2 436.5417
3 394.2917
4 356.3167
5 325.075
6 298.3083
7 283.75
8 274.6083
9 271.0417
10 274.8417
11 284.1167
12 295.15
13 318.4083
14 365.8083
15 398.0833
16 410.4833
17 417.7083
18 426.1
19 424.0583
20 399.8417
21 382.8917
22 385.8583
23 423.45
24 451.7667
/
Wind(t) Scenarios of wind power production
/
1  2
2  1.8583
3  3.1333
4  4.5917
5  4.6
6  16.6417
7  33.475
8  19.2083
9  18.7667
10  48.9417
11  61.0833
12  47.025
13  29.5333
14  28.675
15  55.0417
16  65.825
17  86.15
18  91.8
19  90.8
20  77.55
21  67.8
22  63.65
23  56.8
24  48.1667
/
;
Parameter
D(t) Net power demand (MW);
D(t)=Demand(t)-Wind(t)
;
