# Procedure

1. Run `Inverted_Pendulum9_data.m`
2. Run `parameter_estimation.m`
3. Open `parameter_estimation_2019a.slx`
4. Build and run.
5. Results are saved in `estimation.mat`
6. Stations are ordered in counter-clockwise order in RUCH 337 with the axis pointing upward.

* The order of variables are:
  1. xtilde
  2. xtildedot
  3. mhat
  4. bhat
  5. u
  6. f --> adaptation law for mhat
  7. g --> adaptation law for bhat
