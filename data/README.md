# Demand and wind series

Four CSV files, one row per day and one column per hour, in MW:

| File | Period | Days |
|---|---|---|
| `demand_2019-2021.csv`, `wind_2019-2021.csv` | 19 Dec 2019 20:00 – 19 Dec 2021 20:00 | 731 |
| `demand_2021-2023.csv`, `wind_2021-2023.csv` | 19 Dec 2021 20:00 – 19 Dec 2023 20:00 | 730 |

The first period is used to generate the labels, train the network and calibrate the
fixing bands; the second is the evaluation period of the paper. The net demand, which is
the only input of the network, is the difference between the two, and the day of the
model runs from 20:00 to 20:00.

## Source

Total demand and wind generation of the Gran Canaria system, published by **Red
Eléctrica de España** at <https://demanda.ree.es/visiona/canarias/gcanariaau/total> with
a resolution of 5 minutes. The series distributed here are the hourly averages of those
records, and are redistributed for reproducibility; their terms of use are those of the
source.

## Rebuilding them

`rebuild_profiles.m` regenerates the four files from a `.mat` holding the 5-minute
records in a struct `Data` with fields `ts` (timestamps), `dem` (demand) and `eol`
(wind). Edit the path at the top of the script and run it; the files it writes should be
identical to the ones distributed, to the fourth decimal.
