# Neutron Star Code — Radial Dimensionality Theory Paper III (V8)

Code for **Radial Dimensionality Theory III: Neutron Star Structure under
Solar-Constrained Density-Dependent Dimensionality**
(C. K. Merrill and Claude, 2026). V8 supersedes V7 (its force-rescaling
scripts are archived in [`archive_v7/`](archive_v7/); the earlier V6 code is
in [`archive_v6/`](archive_v6/)).

**What changed in V8 — the sign.** V7 implemented RDT as a force-rescaling
factor F(ρ) = (d_eff − 1)/2 multiplying the TOV pressure gradient, which
weakens gravity and *expands* the star (ΔM_max = +0.9%). That is **not** the
mechanism of Papers I–II, which modify the gravitational **divergence
operator**. The faithful relativistic analog of that operator (this repo)
**contracts** the star: the radial field obeys a (3+δ)-dimensional Gauss law
with matter and the mass integral kept 3D. On the same ~0.1%-validated GR
baseline (M_max = 2.050 M⊙, R(1.4) = 11.69 km), the maximum solar-allowed
effect (A < 0.0087 from Paper I) is a **contraction**:
ΔM_max = −0.71% (Emb1; −0.57% Emb2), ΔR(1.4) = −0.42%, Λ(1.4) reduced. RDT
makes neutron stars *smaller* and *lowers* M_max — matching the Sun/white-dwarf
deficit of Papers I–II, and reversing V7's heavy-pulsar conclusion (RDT does
not help; it slightly worsens the soft-EOS tension).

## Requirements

Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib`.

## Files (V8 — operator form)

| File | Purpose |
|---|---|
| `eos.py` | SLy EOS loaders: CompOSE RG(SLY4) full crust-to-core table (primary) and the Douchin & Haensel (2001) 77-point table (cross-check). Both ε conventions kept — `rho_c2` (correct: tabulated ρ IS the mass-energy density) and `table` (the V6 "recalculated" column, retained as a regression case). |
| `tov.py` | Baseline GR TOV integrator + physical constants used by `tov_operator.py`. (It also still contains the **superseded** V7 F-multiplier `F_rdt`/RDT branch; V8 does not use it — only the GR `solve_star` and constants.) |
| `tov_operator.py` | **The V8 mechanism.** Divergence-operator modified TOV via two inequivalent relativistic embeddings of the (3+δ)-D Gauss law: Embedding 1 (operator correction factor R_op on the acceleration) and Embedding 2 (effective-dimension mass). Both reduce to exact GR TOV at δ=0; A=0 routes to the GR baseline. |
| `validate.py` | GR-baseline gates. **Run first.** (1) Solver vs CompOSE reference (~0.1%); (2) ε-convention regression (correct convention → 2.050; V6 column → 1.793); (3) numerics stability. |
| `gate_independence.py` | **Operator-form gates (the V8 sign-trust check).** Forces the modified path at vanishing A and shows (a) it reproduces GR to integrator level and (b) (ΔM/M)/A is constant across A=1e-5..1e-2 — i.e. the contraction is linear in A and the path → GR as A→0, produced by the same code path with only A varied. (Not the δ=0 fallback.) |
| `boundary_sweep.py` | RDT contraction along Paper I's A₉₅(α) exclusion boundary, both embeddings; prints the ceiling (\|ΔM_max\| ≤ 0.71%, \|ΔR(1.4)\| ≤ 0.42% at the global ceiling α≈0.122, A=0.0087). ~1 min. |
| `figures_v8.py` | Regenerates the paper figure (`figures/fig_mr_v8.pdf`): M–R with the RDT curve *inside* GR (contraction), and the fixed-mass contraction envelope. |
| `data/` | EOS tables: CompOSE RG(SLY4) (`eos.thermo`, `eos.nb`, `eos.mr` reference curve; compose.obspm.fr), D&H 2001 (`eos_sly4_dh2001_v2.dat`), APR (`eos_apr_pc2017.dat`, kept for future multi-EOS work). |
| `derivation.md`, `SIGN_VERDICT.md` | The operator-form derivation (the two embeddings, gates) and the committed sign verdict (audit branch 1a: contraction; numbers, robustness, embedding non-uniqueness). |
| `solar_model.py`, `paper6_alpha1_bound.py` | **Paper VI support** (unchanged): reruns the Paper I solar analysis at the completion's predicted exponent α=1 → A₉₅(α=1) = 0.0120, ε_sat < 0.036, d_eff > 2.964. `solar_model.py` is a verbatim copy of the Paper I V3 model module. ~30 s. |
| `archive_v7/` | Superseded V7 force-rescaling scripts (`predict.py`, `figures.py`, `results_predictions.csv`) — the F-multiplier that carried the wrong sign (+0.9% expansion). Retained for the record. |
| `archive_v6/` | Superseded V6 code/results (the energy-density reconstruction bug documented in V7). |

## Usage

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

python3 validate.py            # GR-baseline gates — run first
python3 gate_independence.py   # operator-form gates: modified path -> GR linearly in A
python3 boundary_sweep.py      # contraction ceiling along A95(alpha)  (~1 min)
python3 figures_v8.py          # paper figure -> figures/fig_mr_v8.pdf
python3 paper6_alpha1_bound.py # Paper VI: solar bound at the predicted exponent (~30 s)
```

This repository serves Papers III and VI of the series: the operator-form
TOV pipeline for Paper III (V8), and the α = 1 solar-bound computation quoted
in Paper VI (which reuses the Paper I machinery via `solar_model.py`).

## Note on the sign reversal

V8 reverses a published headline: V7 reported RDT *raises* M_max and eases
heavy-pulsar tension; the canonical operator form shows RDT *lowers* M_max and
does not help. The reversal is the third self-correction in this paper's line
(after the V1–V2 radius bug and the V6→V7 energy-density fix). The two
embeddings agree in sign and bracket the magnitude at ~0.15 pp — the residual
non-uniqueness of the relativistic embedding of the operator (a magnitude, no
longer a sign, ambiguity), documented in `SIGN_VERDICT.md`. The tidal
deformability reduction is computed in Paper IV's pipeline, not scaled here.

## License

MIT — see [LICENSE](LICENSE).
