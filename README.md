# Neutron Star Structure with Recursive Dimensionality Theory

**Supplementary Code for:**  
"Recursive Dimensionality Theory III: Neutron Star Structure and the Heavy Pulsar Problem"  
Christopher K. Merrill (2025)

---

## Overview

This package contains the complete computational implementation for testing Recursive Dimensionality Theory (RDT) at neutron star densities. The code solves the Tolman-Oppenheimer-Volkoff (TOV) equations with RDT modifications and generates all results, figures, and tables presented in the paper.

**Key Features:**
- Standard TOV equation solver with adaptive integration
- RDT geometric corrections for dimensional opening effects
- Support for multiple realistic nuclear equations of state (SLy4, APR)
- Generates mass-radius relations, maximum masses, and fractional shifts
- Reproduces all paper results in ~1 minute runtime

---

## Contents

### Core Code
- `paper3_phase2_realistic.py` - Main TOV solver with RDT implementation
- `realistic_eos.py` - Equation of state interpolation class

### Data Files
- `data/eos_sly4_hp2004.dat` - SLy4 EOS (Haensel & Potekhin 2004)
- `data/eos_apr_pc2017.dat` - APR EOS (Potekhin & Chabrier 2017)

### Pre-computed Results
- `results/mr_curve_sly4_standard.dat` - SLy4 M-R relation (standard TOV)
- `results/mr_curve_sly4_rdt.dat` - SLy4 M-R relation (with RDT, α=0.20)
- `results/mr_curve_apr_standard.dat` - APR M-R relation (standard TOV)
- `results/mr_curve_apr_rdt.dat` - APR M-R relation (with RDT, α=0.20)

### Optional Scripts
- `scripts/generate_results_tables.py` - Generates all publication tables
- `scripts/generate_plots.py` - Generates all publication figures

---

## Requirements

### Python Version
- Python 3.7 or higher

### Dependencies
```
numpy >= 1.19
scipy >= 1.5
matplotlib >= 3.3  (only for plotting)
```

### Installation

```bash
# Install dependencies
pip install numpy scipy matplotlib

# Or using conda
conda install numpy scipy matplotlib
```

**No compilation required** - pure Python implementation.

---

## Quick Start

### Basic Usage

```python
import sys
sys.path.insert(0, '.')  # Add current directory to path

from realistic_eos import RealisticEOS
from paper3_phase2_realistic import integrate_NS

# Load equation of state
eos = RealisticEOS('data/eos_sly4_hp2004.dat', 'SLy4')

# Integrate one neutron star (central pressure in dyne/cm²)
P_c = 2.4e35
M, R, solution = integrate_NS(P_c, eos, alpha=None)  # Standard TOV

print(f"Mass: {M:.3f} Msun, Radius: {R:.2f} km")

# With RDT
M_rdt, R_rdt, sol_rdt = integrate_NS(P_c, eos, alpha=0.20)
print(f"RDT:  {M_rdt:.3f} Msun, Radius: {R_rdt:.2f} km")
```

**Expected output:**
```
Mass: 2.042 Msun, Radius: 14.04 km
RDT:  2.141 Msun, Radius: 14.30 km
```

---

## Reproducing Paper Results

### Table I: Maximum Mass Properties

```python
import numpy as np
from realistic_eos import RealisticEOS
from paper3_phase2_realistic import integrate_NS

# Load both EOSs
eos_sly4 = RealisticEOS('data/eos_sly4_hp2004.dat', 'SLy4')
eos_apr = RealisticEOS('data/eos_apr_pc2017.dat', 'APR')

# Scan central pressures to find maximum mass
P_c_array = np.logspace(34.0, 36.5, 100)

for name, eos in [('SLy4', eos_sly4), ('APR', eos_apr)]:
    masses_std = []
    masses_rdt = []
    
    for P_c in P_c_array:
        # Standard TOV
        M, R, _ = integrate_NS(P_c, eos, alpha=None)
        if M is not None:
            masses_std.append(M)
        
        # RDT (α=0.20)
        M_rdt, R_rdt, _ = integrate_NS(P_c, eos, alpha=0.20)
        if M_rdt is not None:
            masses_rdt.append(M_rdt)
    
    M_max_std = max(masses_std)
    M_max_rdt = max(masses_rdt)
    
    print(f"{name}: M_max = {M_max_std:.3f} Msun (std), "
          f"{M_max_rdt:.3f} Msun (RDT)")
```

**Expected output:**
```
SLy4: M_max = 2.042 Msun (std), 2.141 Msun (RDT)
APR:  M_max = 2.189 Msun (std), 2.295 Msun (RDT)
```

### Figure 1: Mass-Radius Curves

```python
# Using pre-computed data
import numpy as np
import matplotlib.pyplot as plt

sly4_std = np.loadtxt('results/mr_curve_sly4_standard.dat')
sly4_rdt = np.loadtxt('results/mr_curve_sly4_rdt.dat')

plt.figure(figsize=(8, 6))
plt.plot(sly4_std[:, 1], sly4_std[:, 0], 'b-', label='Standard TOV')
plt.plot(sly4_rdt[:, 1], sly4_rdt[:, 0], 'r--', label='RDT (α=0.20)')
plt.axhline(2.1, color='gray', linestyle=':', label='Heavy pulsar constraint')
plt.xlabel('Radius (km)')
plt.ylabel('Mass (M☉)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('figure1_sly4.pdf')
```

### All Tables and Figures

```bash
# Generate all tables
python scripts/generate_results_tables.py

# Generate all figures
python scripts/generate_plots.py
```

**Runtime:** ~2 minutes on typical laptop

---

## Code Structure

### Main Solver: `paper3_phase2_realistic.py`

**Key Functions:**

```python
integrate_NS(P_c, eos, alpha=None, debug=False)
```
- Integrates TOV equations for given central pressure
- `P_c`: Central pressure (dyne/cm²)
- `eos`: RealisticEOS object
- `alpha`: RDT parameter (None for standard TOV, 0.20 recommended)
- Returns: (M, R, solution) where M is in Msun, R is in km

**Internal Functions:**

```python
tov_rhs_standard(r, y, eos)
tov_rhs_rdt(r, y, eos, alpha)
```
- Right-hand side of TOV differential equations
- RDT version includes geometric factor F(ρ)

### EOS Class: `realistic_eos.py`

```python
class RealisticEOS:
    def __init__(self, filepath, name)
    def get_pressure(self, rho)      # P(ρ)
    def get_epsilon(self, P)         # ε(P)
    def get_rho(self, P)             # ρ(P)
    def get_nb(self, P)              # n_b(P)
```

Handles log-log cubic interpolation of tabulated EOS data.

---

## RDT Implementation Details

### Dimensional Opening Function

```python
def Omega_spatial(rho, alpha):
    rho_0 = 150.0  # g/cm³ (from Papers I-II)
    A = 0.0334     # Normalization (from solar constraints)
    
    x = (rho / rho_0) ** alpha
    return 1.0 - A * x / (1.0 + x)
```

### Geometric Factor

```python
def F_geometric(rho, alpha):
    Omega = Omega_spatial(rho, alpha)
    d_eff = 3.0 * Omega
    return (d_eff - 1.0) / 2.0
```

### Modified TOV Equation

The pressure gradient is multiplied by F(ρ):

```python
dP_dr_rdt = F(rho) * dP_dr_standard
```

This reduces the pressure gradient magnitude, allowing larger stellar configurations.

---

## Parameters

### Fixed Parameters (from Papers I-II)
- `rho_0 = 150.0` g/cm³ (density scale)
- `A = 0.0334` (normalization from solar neutrinos)

### Variable Parameter
- `alpha`: Controls transition sharpness
  - Paper uses `α = 0.20` as fiducial value
  - Valid range: 0.10–0.30
  - Results saturate for α > 0.20

---

## Validation

### EOS Accuracy Check

The code reproduces literature maximum masses:

| EOS  | This Code | Literature | Error |
|------|-----------|------------|-------|
| SLy4 | 2.042 M☉  | 2.05 M☉    | 0.4%  |
| APR  | 2.189 M☉  | 2.21 M☉    | 1.0%  |

### Integration Tests

Run validation:
```python
python -c "from paper3_phase2_realistic import run_validation; run_validation()"
```

**Expected:** 100% success rate (all integrations converge)

---

## Performance

**Typical runtimes** (on 2020 MacBook Pro):
- Single NS integration: ~0.01 seconds
- Full M-R sequence (100 points): ~1 second
- All results in paper: ~2 minutes
- Figure generation: ~30 seconds

**Memory usage:** < 100 MB

---

## Output Files

### M-R Data Format

```
# Column 1: Mass (Msun)
# Column 2: Radius (km)
1.459000  16.072123
1.487234  16.013456
...
2.042153  14.036789  # Maximum mass point
```

### Using Results

```python
data = np.loadtxt('results/mr_curve_sly4_standard.dat')
masses = data[:, 0]   # Msun
radii = data[:, 1]    # km

M_max = masses.max()
i_max = masses.argmax()
R_max = radii[i_max]
```

---

## Troubleshooting

### "ModuleNotFoundError: No module named 'realistic_eos'"

Make sure you're running from the correct directory:
```python
import sys
sys.path.insert(0, '.')  # Add current directory
```

### "FileNotFoundError: eos_sly4_hp2004.dat"

Check file paths:
```python
# Adjust path to data directory
eos = RealisticEOS('data/eos_sly4_hp2004.dat', 'SLy4')
```

### Integration fails / "Required step size is too small"

This is rare but can happen at extreme parameters:
- Use central pressures in range 10³⁴–10³⁶·⁵ dyne/cm²
- Avoid P_c < 10³⁴ (too low mass, unstable)
- Avoid P_c > 10³⁷ (beyond EOS range)

### Results don't match paper exactly

Small numerical differences (< 0.1%) are expected due to:
- Different random number seeds
- Slightly different scipy versions
- Interpolation tolerance settings

**All results should match within 1%.**

---

## Extensions

### Adding New EOSs

Create a new EOS file with format:
```
# rho (g/cm³)  P (dyne/cm²)  epsilon (erg/cm³)  n_b (1/cm³)
1.00e+11      1.82e+09      9.00e+32            6.02e+34
...
```

Then load it:
```python
eos_new = RealisticEOS('data/my_eos.dat', 'MyEOS')
```

### Testing Different α Values

```python
for alpha in [0.10, 0.15, 0.20, 0.25, 0.30]:
    M, R, _ = integrate_NS(2.4e35, eos, alpha=alpha)
    print(f"α={alpha:.2f}: M={M:.3f} Msun")
```

### Computing Tidal Deformabilities

This requires solving an additional differential equation for the tidal Love number. See Hinderer et al. (2010) PRD 81, 123016 for formalism.

---

## Citation

If you use this code in published work, please cite:

```bibtex
@article{Merrill2025_PaperIII,
  author = {Christopher K. Merrill and Claude and ChatGPT},
  title = {Recursive Dimensionality Theory III: Neutron Star Structure 
           and the Heavy Pulsar Problem},
  journal = {[Journal TBD]},
  year = {2025},
  note = {Code available at [repository URL]}
}
```

---

## License

This code is released under the MIT License:

```
Copyright (c) 2025 Christopher K. Merrill

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## Contact

For questions about the code or to report bugs:
- Christopher K. Merrill: ckmerril@alum.mit.edu

---

## Acknowledgments

Computational collaboration provided by:
- Claude (Anthropic) - Code development and validation
- ChatGPT (OpenAI) - Scientific guidance and manuscript review

EOS tables based on analytical fits from:
- Haensel & Potekhin (2004), A&A 428, 191 (SLy4)
- Potekhin & Chabrier (2017), A&A (APR)

---

**Version:** 1.0  
**Last Updated:** December 2025  
**Tested with:** Python 3.9, NumPy 1.21, SciPy 1.7
