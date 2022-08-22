# crossing_statistics
Final project for the 2022 Serrapilheira QBio Training Program Computation Methods Course.

In this project, we efficiently generate trajectories of an two-dimensional Ornstein-Uhlenbeck process and compile their first-encounter statistics. To do so, we implement two Julia modules for pseudo-random generation and slightly generalizable integration for linear stochastic differentiable equations. We then generate a dataset of encounter events and visualize our data using a histogram of the first-encounter times and an animation of the time-evolution of the motion.

## Project structure

```
project/
*    ├── data/
     ├── docs/
*    ├── figs/
     ├── scripts/
*    └── README.md
```

The `scripts/` folder contains three scripts to be run using Julia 1.8, as well as the modules `PRNG.jl`, for pseudo-random number generation, and `SDE.jl` for linear stochastic differential equation integration. To run the scripts, one simply needs to call

`julia ./scripts/script.jl`

In the `docs/` folder one finds the `report.Rmd` R Markdown file, compiled into a `.pdf` using `knitr` with `pdflatex`.

## Required packages

The Julia scripts require base Julia >=1.8, as well as the `DelimitedFiles` and `Plots` package. The `JuliaCall` package on `R` is required for compiling the report. 