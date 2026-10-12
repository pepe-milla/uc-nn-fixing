* Variables

Positive variables

*Variables regarding the TGUs
p(i,q,t)                 Output power of TGU i corresponding to the qth segment of the production cost curve [MW]
pv(i,t)                  Power variation between two consecutive hours [MW]

*Variables regarding the PHS
pg(c,t)                  Output power over Pmin of unit c in period t (turbine mode) [MW]
pp(c,t)                  Input power over Pmin of unit c in period t (pump mode) [MW]
pvp(c,t)                 Power variation between two consecutive hours [MW]

*General variables
wc(t)                    Wind curtailment in period t [MW]
v(t)                     volume in the upper reservoir in period t [m3]
;

*Upper limits of some varaiables
wc.up(t)=Wind(t);
v.up(t)=3770000;

Free variables
z                        Objective function - scheduled power supply costs(�)

Binary variables

*Variables regarding the TGUs
u(i,t)                   Commitment of the TGU i during period t
delta(i,t)               Binary variable which takes the value 1 if the TGU i starts-up in period t and 0 otherwise
w(i,t)                   Binary variable which takes the value 1 if the TGU i shuts-down in period t and 0 otherwise

*Variables regarding the PHS
ug(c,t)                  Commitment of the PSH unit c in turbine mode during period t
deltag(c,t)              Binary variable which takes the value 1 if the PSH unit c starts-up in turbine mode in period t and 0 otherwise
up(c,t)                  Commitment of the PSH unit c in pump mode during period t
deltap(c,t)              Binary variable which takes the value 1 if the PSH unit c starts-up in pump mode in period t and 0 otherwise
;

