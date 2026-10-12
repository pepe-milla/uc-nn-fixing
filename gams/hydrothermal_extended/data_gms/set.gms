* Data reading from Matlab
options LIMROW=7500, LIMCOL=1500
Set
*set timing data
t   hours to consider in the calculation /1*24/
i   thermal generating units /1*16/
c   pum-turbine units /1*6/
q segments of the hourly fuel cost curve of the thermal units /1*10/
;
alias(t,tt);
alias(i,ii);
alias(c,cc);
