# Neural-network-assisted MIP for the Unit Commitment problem

Models, code, trained networks and results of the paper

> J. Milla, J. I. Pérez-Díaz, D. Fernández, *Application of neural networks to assist
> MIP solvers in the UC in power systems with penetration of renewable generation and
> pumped-storage*.

The case study is the isolated power system of Gran Canaria (Spain): 16 thermal units
and the pumped-storage hydropower plant (PSHP) under construction. A neural network
predicts the commitment of the units from the net demand profile of the day, and the
predictions it is confident about are fixed before the MILP is handed to the solver.

Four test cases cross the two systems with the two formulations:

| Case | System | Formulation | Stopping criteria |
|---|---|---|---|
| `thermal_base` | thermal only | base | 0.1 %, 120 s |
| `thermal_extended` | thermal only | + minimum up/down times and ramps | 0.1 %, 120 s |
| `hydrothermal_base` | with PSHP | base | 1.5 %, 1200 s |
| `hydrothermal_extended` | with PSHP | + minimum up/down times and ramps | 1.5 %, 1200 s |

## Reproducing the tables and figures

Everything reported in the paper can be regenerated from the stored results, without
solving a single MILP:

```matlab
cd reproduce
reproduce_all
```

It takes a couple of minutes, prints Tables 5 to 14 and writes Figures 3 and 4 into
`reproduce/figures`. It needs MATLAB (R2023b or later, because the networks are
`dlnetwork` objects trained with `trainnet`) with the Deep Learning and Statistics
toolboxes. No GAMS licence is required for this.

## Re-running the campaigns

This does need GAMS and CPLEX, which are commercial and are not distributed here. The
campaigns of the paper were run with **GAMS 23.5 and CPLEX 12.2**, in deterministic
parallel mode with two threads per solver process. Before running anything, edit
`matlab/common/uc_config.m`, which is the only file with machine-dependent paths:

```matlab
c.gamsExe    = 'C:\Program Files\GAMS23.5\gams.exe';
c.workerRoot = fullfile(tempdir, 'uc_workers');
```

The order for one case is: generate the labels, train the network, calibrate the bands,
run the unassisted baseline, run the assisted campaign. For example, for the
hydrothermal extended case:

```matlab
cd matlab/hydrothermal_extended/scripts
generate_training_data_par      % 731 days of 2019-2021 -> results/MIP_Results_2019-2021.mat
define_NN                       % -> trained models/NN.mat
calibrate_fixing                % -> data/Fixing_bands.mat
generate_baseline_par           % 730 days of 2021-2023, unassisted
assess_NN_par                   % the same days with the fixings of Strategy 2
```

Be warned about the scale: the four baselines and the assisted campaigns of the paper
took about 220 hours of solver time on a workstation with an Intel Core i9-13900K,
solving four days at a time.

In `hydrothermal_base` there is no separate baseline script: the unassisted baseline was
produced with `generate_training_data_par.m`, changing the period, the time limit and
the name of the output file, as its header explains. The labels of that case were
generated with the same script and a time limit of 600 s; every campaign reported in the
paper uses 1200 s.

The solver statistics of Table 12 are obtained by re-solving with the CPLEX log kept:
`matlab/solver_metrics/scripts/run_solver_metrics.m`, summarised by `sm_summary.m`.

## Layout

```
data/                    hourly demand and wind series of both periods (CSV)
gams/<case>/             the GAMS model: run file, data_gms/, cplex.op9
matlab/common/           uc_config.m (paths) and uc_load_profiles.m (series)
matlab/<case>/scripts/   label generation, training, calibration and campaigns
matlab/<case>/data/      labels and calibrated bands
matlab/<case>/results/   campaign results (.mat)
matlab/<case>/trained models/   the networks
matlab/solver_metrics/   re-solves keeping the CPLEX log, and parses it
reproduce/               regenerates every table and figure of the paper
```

The `.mat` files of `results/` hold, for every day: the commitment and output of every
unit, the cost, the optimality gap, the computation time and, for the assisted runs,
the number of fixed variables.

## Notes

- **Curtailment penalty.** The base formulation penalises wind curtailment at
  10<sup>6</sup> €/MWh, which in practice forbids it; the extended formulation sets the
  penalty to zero, because the minimum up and down times make curtailment unavoidable on
  some days. Costs are therefore comparable within a formulation but not across them.
- **Initial conditions.** Every day starts from the same initial status and output of
  the units, and the same reservoir volume, so that all days and cases are comparable.
  In the thermal extended case this makes days 660 and 661 of the evaluation period
  infeasible; they are excluded from the results, as explained in the paper.
- **Computation times.** In the thermal extended case the times stored by the first
  campaigns were re-read from the solver logs with a corrected parser
  (`matlab/thermal_extended/scripts/fix_comp_time.m`); the files distributed here
  contain the corrected values.

## Data

The demand and wind series come from Red Eléctrica de España. See `data/README.md`.

## Licence

Code and models: MIT (see `LICENSE`), copyright Universidad Politécnica de Madrid,
by José Milla López-Asiaín, Juan I. Pérez-Díaz and Daniel Fernández-Muñoz. The
parameters of the Gran Canaria system are those published in the paper and in the
reference it takes them from; the data of Red Eléctrica de España keep their own
terms of use.
