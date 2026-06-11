"""
Generate Results Tables for Paper III

Create publication-ready tables showing:
1. Standard TOV results
2. RDT effects
3. Key radius values at M=1.4 Msun
"""

import sys
sys.path.insert(0, '/home/claude')
sys.path.insert(0, '/mnt/user-data/outputs')
from realistic_eos import RealisticEOS
from paper3_phase2_realistic import integrate_NS
import numpy as np

print("="*70)
print("PAPER III RESULTS TABLES")
print("="*70)

# Load EOS
eos_sly4 = RealisticEOS('/mnt/user-data/outputs/eos_sly4_hp2004.dat', 'SLy4')
eos_apr = RealisticEOS('/mnt/user-data/outputs/eos_apr_pc2017.dat', 'APR')

# Alpha for RDT (use 0.20 as representative value per ChatGPT suggestion)
alpha_rep = 0.20

# Generate M-R sequences
P_c_array = np.logspace(34.0, 36.5, 150)  # Extended range, higher resolution

results = {}

for eos_name, eos in [('SLy4', eos_sly4), ('APR', eos_apr)]:
    print(f"\nProcessing {eos_name}...")
    
    results[eos_name] = {'Standard': {}, 'RDT': {}}
    
    for alpha_label, alpha_val in [('Standard', None), ('RDT', alpha_rep)]:
        masses = []
        radii = []
        P_centers = []
        
        for P_c in P_c_array:
            M, R, sol = integrate_NS(P_c, eos, alpha=alpha_val, debug=False)
            if M is not None and M > 0.5:
                masses.append(M)
                radii.append(R)
                P_centers.append(P_c)
        
        masses = np.array(masses)
        radii = np.array(radii)
        P_centers = np.array(P_centers)
        
        # Find maximum mass
        i_max = np.argmax(masses)
        M_max = masses[i_max]
        R_max = radii[i_max]
        P_c_max = P_centers[i_max]
        
        # Find R(1.4 Msun) - interpolate
        if masses.min() < 1.4 < masses.max():
            R_1p4 = np.interp(1.4, masses, radii)
        else:
            R_1p4 = None
        
        results[eos_name][alpha_label] = {
            'masses': masses,
            'radii': radii,
            'P_centers': P_centers,
            'M_max': M_max,
            'R_max': R_max,
            'R_1p4': R_1p4,
            'P_c_max': P_c_max
        }

# Print Table 1: Maximum Mass Properties
print("\n" + "="*70)
print("TABLE 1: MAXIMUM MASS PROPERTIES")
print("="*70)
print("\n| EOS   | Model    | M_max (Msun) | R_max (km) | P_c (dyne/cm²) |")
print("|-------|----------|--------------|------------|----------------|")

for eos_name in ['SLy4', 'APR']:
    # Standard
    M = results[eos_name]['Standard']['M_max']
    R = results[eos_name]['Standard']['R_max']
    P = results[eos_name]['Standard']['P_c_max']
    print(f"| {eos_name:5} | Standard | {M:12.3f} | {R:10.2f} | {P:14.3e} |")
    
    # RDT
    M_rdt = results[eos_name]['RDT']['M_max']
    R_rdt = results[eos_name]['RDT']['R_max']
    P_rdt = results[eos_name]['RDT']['P_c_max']
    print(f"| {eos_name:5} | RDT α=0.2 | {M_rdt:12.3f} | {R_rdt:10.2f} | {P_rdt:14.3e} |")

# Print Table 2: R(1.4 Msun) 
print("\n" + "="*70)
print("TABLE 2: CANONICAL 1.4 Msun NEUTRON STAR")
print("="*70)
print("\n| EOS   | Model    | M (Msun) | R (km) |")
print("|-------|----------|----------|--------|")

for eos_name in ['SLy4', 'APR']:
    R_std = results[eos_name]['Standard']['R_1p4']
    R_rdt = results[eos_name]['RDT']['R_1p4']
    
    if R_std is not None:
        print(f"| {eos_name:5} | Standard | {1.4:8.2f} | {R_std:6.2f} |")
    if R_rdt is not None:
        print(f"| {eos_name:5} | RDT α=0.2 | {1.4:8.2f} | {R_rdt:6.2f} |")

# Print Table 3: RDT Effects (Fractional Changes)
print("\n" + "="*70)
print("TABLE 3: RDT EFFECTS (α = 0.20)")
print("="*70)
print("\n| EOS   | ΔM_max (Msun) | ΔM_max (%) | ΔR_max (km) | ΔR_max (%) | ΔR(1.4) (km) | ΔR(1.4) (%) |")
print("|-------|---------------|------------|-------------|------------|--------------|-------------|")

for eos_name in ['SLy4', 'APR']:
    M_std = results[eos_name]['Standard']['M_max']
    M_rdt = results[eos_name]['RDT']['M_max']
    R_std = results[eos_name]['Standard']['R_max']
    R_rdt = results[eos_name]['RDT']['R_max']
    
    dM = M_rdt - M_std
    dM_pct = 100 * dM / M_std
    dR = R_rdt - R_std
    dR_pct = 100 * dR / R_std
    
    R_1p4_std = results[eos_name]['Standard']['R_1p4']
    R_1p4_rdt = results[eos_name]['RDT']['R_1p4']
    
    if R_1p4_std is not None and R_1p4_rdt is not None:
        dR_1p4 = R_1p4_rdt - R_1p4_std
        dR_1p4_pct = 100 * dR_1p4 / R_1p4_std
    else:
        dR_1p4 = 0
        dR_1p4_pct = 0
    
    print(f"| {eos_name:5} | {dM:13.3f} | {dM_pct:10.2f} | {dR:11.2f} | {dR_pct:10.2f} | {dR_1p4:12.2f} | {dR_1p4_pct:11.2f} |")

# Save M-R data for plotting
print("\n" + "="*70)
print("SAVING M-R DATA FOR PLOTTING")
print("="*70)

for eos_name in ['SLy4', 'APR']:
    for model in ['Standard', 'RDT']:
        masses = results[eos_name][model]['masses']
        radii = results[eos_name][model]['radii']
        
        data = np.column_stack([masses, radii])
        filename = f'/mnt/user-data/outputs/mr_curve_{eos_name.lower()}_{model.lower()}.dat'
        
        header = f"""M-R curve for {eos_name} ({model})
Columns: Mass (Msun), Radius (km)
Model: {'Standard TOV' if model == 'Standard' else f'RDT with α={alpha_rep}'}
Generated: 2025-12-05
"""
        np.savetxt(filename, data, fmt='%.6f', header=header)
        print(f"Saved: {filename}")

# Summary statistics
print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

# EOS universality check
dM_sly4 = 100 * (results['SLy4']['RDT']['M_max'] - results['SLy4']['Standard']['M_max']) / results['SLy4']['Standard']['M_max']
dM_apr = 100 * (results['APR']['RDT']['M_max'] - results['APR']['Standard']['M_max']) / results['APR']['Standard']['M_max']
dR_sly4 = 100 * (results['SLy4']['RDT']['R_max'] - results['SLy4']['Standard']['R_max']) / results['SLy4']['Standard']['R_max']
dR_apr = 100 * (results['APR']['RDT']['R_max'] - results['APR']['Standard']['R_max']) / results['APR']['Standard']['R_max']

print(f"\nFractional mass shifts:")
print(f"  SLy4: {dM_sly4:.2f}%")
print(f"  APR:  {dM_apr:.2f}%")
print(f"  Difference: {abs(dM_sly4 - dM_apr):.2f}% (EOS universality)")

print(f"\nFractional radius shifts:")
print(f"  SLy4: {dR_sly4:.2f}%")
print(f"  APR:  {dR_apr:.2f}%")
print(f"  Difference: {abs(dR_sly4 - dR_apr):.2f}% (EOS universality)")

# Heavy pulsar test
print(f"\nHeavy pulsar constraint (M_max > 2.1 Msun):")
print(f"  SLy4 Standard: {results['SLy4']['Standard']['M_max']:.3f} Msun {'✗ FAILS' if results['SLy4']['Standard']['M_max'] < 2.1 else '✓ PASSES'}")
print(f"  SLy4 + RDT:    {results['SLy4']['RDT']['M_max']:.3f} Msun {'✗ FAILS' if results['SLy4']['RDT']['M_max'] < 2.1 else '✓ PASSES'}")
print(f"  APR Standard:  {results['APR']['Standard']['M_max']:.3f} Msun {'✗ FAILS' if results['APR']['Standard']['M_max'] < 2.1 else '✓ PASSES'}")
print(f"  APR + RDT:     {results['APR']['RDT']['M_max']:.3f} Msun {'✗ FAILS' if results['APR']['RDT']['M_max'] < 2.1 else '✓ PASSES'}")

print("\n" + "="*70)
print("RESULTS TABLES COMPLETE")
print("="*70)
print("\nReady for Paper III manuscript!")
