# Radial Dimensionality Theory III: Neutron Star Structure under Solar-Constrained Density-Dependent Dimensionality

Scientific source and reproducible results for version v9.

- [Paper](reproduction/output/pdf/publication/paper3/merrill_RDT_3_Neutrons_and_Pulsars_v9.pdf)
- [Manuscript](reproduction/paper3_NeutronStars/V8_Paper/paper3_v8.tex)

## Reproduction

Requires Python 3.14, NumPy 2.4.2, SciPy 1.17.1, Matplotlib 3.10.8,
TeX Live with revtex4-2, and Poppler. The resource controller uses POSIX
nice/setpriority, SIGPROF and SIGSTOP/SIGCONT.

From the repository directory:

```sh
nice -n 19 python3 reproduce.py --output /tmp/merrill2025_neutronstar_code-reproduction
```

The output directory must not exist. This verifies every capsule checksum,
regenerates tables, figures or algebraic checks from retained scientific outputs,
and builds a PDF whose extracted text matches the included paper.
No parameter scan runs in this command. Full scientific runs use the producing
scripts and their data under `reproduction/`. Use one worker, verified nice 19,
at most 95% of one CPU and one native-library thread.

`release.json` identifies the current manuscript and PDF;
`reproduction/CAPSULE.json` binds the scientific source and data by checksum.
Shared dependencies and saved results are included so no sibling checkout is
required. Historical top-level scripts and data retain their original identity;
the current paper and source manifest distinguish them from the revised results.

## Scientific scope

This version retains the published title. It uses the full density law and corrected central seeds, comparing two specified hydrostatic prescriptions for one equation of state over a sampled parameter domain. The reported mass extrema occur at the equation-of-state domain boundary and do not establish shifts in the true maximum stable mass. The conditional tidal example is reconciled with Paper IV, with the finite cutoff and exterior-matching limitations explicit.
