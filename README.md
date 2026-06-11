# Neutron Star Code — Radial Dimensionality Theory Paper III (V7)

Code for **Radial Dimensionality Theory III: Neutron Star Structure under
Solar-Constrained Density-Dependent Dimensionality**
(C. K. Merrill and Claude, 2026). Supersedes the V6 code, now in
[`archive_v6/`](archive_v6/).

TOV solver for SLy neutron stars (GR baseline + RDT-modified), validated
to ~0.1% against the published CompOSE mass–radius relation. RDT enters as
an interior pressure-gradient rescaling F(ρ) = (d_eff − 1)/2 with the
series opening law, bounded by Paper I's solar constraint (A < 0.0087,
95% CL). Headline ceilings: ΔM_max ≤ +0.9%, ΔR(1.4) ≤ +0.5%,
ΔΛ(1.4) ≤ +2.6%.

## Requirements

Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib`.

## Files

| File | Purpose |
|---|---|
| `eos.py` | SLy EOS loaders: CompOSE RG(SLY4) full crust-to-core table (primary) and the Douchin & Haensel (2001) 77-point table (cross-check). Both ε conventions kept — `rho_c2` (correct: tabulated ρ IS the mass-energy density) and `table` (the V6 "recalculated" column, retained as a regression case). |
| `tov.py` | TOV integrator (GR and RDT), series initial conditions, surface event detection. |
| `validate.py` | Validation gates. **Run first.** (1) Solver vs CompOSE reference curve (<1%; achieves ~0.1%); (2) ε-convention regression — correct convention gives M_max = 2.050, the V6 column reproduces the V6 bug (1.793); (3) numerics stability. |
| `predict.py` | RDT shifts along Paper I's solar exclusion boundary A₉₅(α). Prints the paper-numbers block; writes `results_predictions.csv`. |
| `figures.py` | Regenerates the paper figure (`figures/fig_mr_v7.pdf`). |
| `data/` | EOS tables: CompOSE RG(SLY4) (`eos.thermo`, `eos.nb`, `eos.mr` reference curve; compose.obspm.fr), D&H 2001 (`eos_sly4_dh2001_v2.dat`), APR (`eos_apr_pc2017.dat`, unused in V7 — kept for future multi-EOS work). |
| `archive_v6/` | Superseded V6 code and results, kept for the record. Its energy-density reconstruction (in `realistic_eos.py`) is the bug documented in the V7 paper; its `results/` are the V6 data files. |

## Usage

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

python3 validate.py   # all gates must PASS first
python3 predict.py    # paper-numbers block + CSV  (~1 min)
python3 figures.py    # paper figure
```

## License

MIT — see [LICENSE](LICENSE).
