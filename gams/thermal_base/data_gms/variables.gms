* Variables

Positive variables

*Variables regarding the TGUs
p(i,q,t)                 Output power of TGU i corresponding to the qth segment of the production cost curve [MW]
pv(i,t)                  Power variation between two consecutive hours [MW]

*General variables
wc(t)                    Wind curtailment in period t [MW]
;

*Upper limits of some varaiables
wc.up(t)=Wind(t);

Free variables
z                        Objective function - scheduled power supply costs(€)

Binary variables

*Variables regarding the TGUs
u(i,t)                   Commitment of the TGU i during period t
delta(i,t)               Binary variable which takes the value 1 if the TGU i starts-up in period t and 0 otherwise
;

