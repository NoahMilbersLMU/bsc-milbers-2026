# Stochastic Spectral Sampling for Radiative Transfer

Source code for the bachelor thesis *Stochastic Spectral Sampling for
Radiative Transfer* (Noah Milbers, Meteorological Institute, LMU Munich, 2026;
supervisor Prof. Dr. Bernhard Mayer).

The model is a 1D radiative-convective equilibrium model in C. It solves the
Schwarzschild equation on 100 001 thermal-infrared wavelengths for CO2, H2O,
CH4, N2O and O3. The variants below replace the full line-by-line (LBL)
spectral integration with K sampled wavelengths per time step.

## Contents

| File (`src/`)                            | Method                                      | Thesis   |
|------------------------------------------|---------------------------------------------|----------|
| `spectralCO2H2O.c`                       | Full LBL reference (no sampling)            | Ch. 2    |
| `spectralCO2H2O_shuffle_v1.c`            | Shuffle-based sampling, fixed draw          | Sec. 3.1 |
| `spectralCO2H2O_shuffle_v2.c`            | Shuffle-based sampling, fresh draw per step | Sec. 3.1 |
| `spectralCO2H2O_montecarlo.c`            | Crude Monte Carlo integration               | Sec. 3.2 |
| `spectralCO2H2O_equidistant.c`           | Rotating equidistant sampling               | Sec. 3.3 |
| `spectralCO2H2O_equidistant_autostop.c`  | Rotating equidistant with stopping criterion| Ch. 5    |
| `spectralCO2H2O_importance.c`            | Planck-weighted importance sampling         | Sec. 3.4 |
| `spectralCO2H2O_stratified_v1.c`         | Stratified sampling, uniform bins           | Sec. 3.5 |
| `spectralCO2H2O_stratified_v2.c`         | Stratified sampling, Planck-weighted bins   | Sec. 3.5 |
| `spectralCO2H2O_adaptive_v1.c`           | Adaptive suppressed sampling, variant 1     | Sec. 3.6 |
| `spectralCO2H2O_adaptive_v2.c`           | Adaptive suppressed sampling, variant 2     | Sec. 3.6 |
| `runlog.h`                               | Run configuration and per-step CSV logging  | -        |

## Files not included - request from Prof. Mayer's group

The model builds on code from Prof. Dr. Bernhard Mayer's group (Meteorological
Institute, LMU Munich), which is not redistributed here:

- `ascii.c`, `ascii.h` - ASCII I/O library from libRadtran
  (Kylling & Mayer, <http://www.libradtran.org>), used to read the
  atmosphere and optical-thickness tables
- `schwarzschild.c`, `schwarzschild.h` - Schwarzschild equation solver

Please request these from Prof. Mayer's group. Place them in `src/` next to the
model sources.

## Input data

The model reads its input from `./lbl.arts/` relative to the working directory:

- `fpda.lbl.atm` - atmosphere profile (levels with pressure, temperature and
  gas volume mixing ratios)
- `lbl.co2.asc`, `lbl.h2o.asc`, `lbl.ch4.asc`, `lbl.n2o.asc`, `lbl.o3.asc` -
  per-gas optical thickness tables. Column 1 is the wavelength grid
  (100 001 points, 4 um - 1 mm, 0.025 cm^-1 spacing in wavenumber), the
  other columns give the optical thickness of each layer.

These tables are not included because they are large (~25-50 MB per gas and
atmosphere). The thesis tables were computed with HAPI, the HITRAN
Application Programming Interface, from the HITRAN2024 line list, plus the
MT_CKD 4.3 water-vapour continuum for H2O:

- HAPI documentation: <https://hitran.org/hapi/>
  (code: <https://github.com/hitranonline/hapi>; Kochanov et al. 2016,
  JQSRT 177, 15-30)
- MT_CKD: <https://github.com/AER-RC/MT_CKD_H2O>

## Build

With the files above in `src/`:

```bash
cd src
gcc -O3 -o spectralCO2H2O_montecarlo spectralCO2H2O_montecarlo.c ascii.c schwarzschild.c -I . -lm
```

The same command works for every variant. Only `math.h` and the C standard
library are needed.

## Run

```bash
ln -s /path/to/tables lbl.arts      # directory holding the .atm and .asc files

./spectralCO2H2O <h2o 0/1> <CO2 multiplier>                         # LBL reference
./spectralCO2H2O_montecarlo <h2o 0/1> <CO2 mult> [K] [seed]
./spectralCO2H2O_equidistant <h2o 0/1> <CO2 mult> [K]
./spectralCO2H2O_equidistant_autostop <h2o 0/1> <CO2 mult> [K] [ncycles] [eps]
./spectralCO2H2O_adaptive_v1 <h2o 0/1> <CO2 mult> [K] [seed] [suppress] [rectime]
```

The shuffle, importance and stratified variants take the same arguments as
`montecarlo`. `K = 0` runs the full LBL spectrum. If no seed is given, the
current time is used. Each binary prints its usage when called without
arguments.

Optional environment variables (defined in `runlog.h`). With none set, the
model runs exactly as described above:

| Variable      | Effect                                                 | Default                 |
|---------------|--------------------------------------------------------|-------------------------|
| `LBL_ATMFILE` | path of the atmosphere profile                         | `lbl.arts/fpda.lbl.atm` |
| `LBL_SIMDAYS` | simulated duration in days                             | `1095` (3 years)        |
| `LBL_OUTFILE` | write a per-step CSV `step,day,dt_s,Tsfc,Ttoa`         | off                     |
| `LBL_META`    | comment header written into the CSV                    | -                       |
