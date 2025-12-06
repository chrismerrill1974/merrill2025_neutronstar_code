"""
Paper III Phase 2: Realistic EOS TOV Solver with RDT Modification

This code implements the Tolman-Oppenheimer-Volkoff (TOV) equations for 
realistic nuclear equations of state (SLy4, APR4) with Recursive Dimensionality 
Theory (RDT) geometric corrections.

Implementation Agreement: Section 3.1.2 - Realistic EOS Implementation
Authors: Christopher K. Merrill, ChatGPT, Claude
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Import the realistic EOS class
sys.path.insert(0, '/home/claude')
sys.path.insert(0, '/mnt/user-data/outputs')
from realistic_eos import RealisticEOS

# ============================================================================
# PHYSICAL CONSTANTS (CGS)
# ============================================================================

G = 6.67430e-8          # gravitational constant [cm^3 g^-1 s^-2]
c = 2.99792458e10       # speed of light [cm/s]
Msun = 1.98847e33       # solar mass [g]
km = 1.0e5              # kilometer [cm]
m_baryon = 1.66e-24     # average baryon mass [g]

# ============================================================================
# RDT PARAMETERS (from Papers I & II)
# ============================================================================

rho_sun_core = 150.0              # solar core density [g/cm^3]
Omega_sun = 2.95 / 3.0            # solar opening fraction
A_saturating = 2.0 * (1.0 - Omega_sun)  # saturation parameter
rho0 = rho_sun_core               # reference density

# ============================================================================
# RDT OPENING FRACTION AND GEOMETRIC FACTOR
# ============================================================================

def Omega_spatial(rho, alpha):
    """
    Saturating dimensional opening fraction.
    
    Parameters:
        rho: baryon density [g/cm^3]
        alpha: density-scaling exponent
    
    Returns:
        Opening fraction Omega in [0, 1]
    """
    x = (rho / rho0)**alpha
    return 1.0 - A_saturating * x / (1.0 + x)

def d_eff_spatial(rho, alpha):
    """
    Effective spatial dimension.
    
    Returns:
        d_eff = 3 * Omega_spatial(rho, alpha)
    """
    return 3.0 * Omega_spatial(rho, alpha)

def F_geom(rho, alpha):
    """
    Geometric correction factor for RDT.
    
    This multiplies the RHS of the TOV equation (Option B).
    
    Returns:
        F = (d_eff - 1) / 2
    """
    d_eff = d_eff_spatial(rho, alpha)
    return (d_eff - 1.0) / 2.0

# ============================================================================
# TOV EQUATIONS WITH REALISTIC EOS
# ============================================================================

def tov_rhs_standard(r, y, eos):
    """
    Standard TOV equations (no RDT) with realistic EOS.
    
    Args:
        r: radial coordinate [cm]
        y: [P, m] = [pressure, enclosed mass]
        eos: RealisticEOS object
    
    Returns:
        [dP/dr, dm/dr]
    """
    P, m = y
    
    # Handle boundary conditions  
    if P <= 0:  # Changed from 1e20 - was stopping integration too early!
        return [0.0, 0.0]
    
    if r <= 0:
        return [0.0, 0.0]
    
    # Get thermodynamic quantities from EOS
    epsilon = eos.get_epsilon(P)
    rho = eos.get_rho(P)
    
    if rho <= 0 or epsilon <= 0:
        return [0.0, 0.0]
    
    # TOV equation
    numerator = G * (epsilon + P) * (m + 4.0*np.pi*r**3 * P/c**2)
    denominator = c**2 * r**2 * (1.0 - 2.0*G*m/(r*c**2))
    
    if denominator <= 0:
        return [0.0, 0.0]
    
    dP_dr = -numerator / denominator
    dm_dr = 4.0*np.pi*r**2 * epsilon / c**2
    
    if not (np.isfinite(dP_dr) and np.isfinite(dm_dr)):
        return [0.0, 0.0]
    
    return [dP_dr, dm_dr]

def tov_rhs_rdt(r, y, eos, alpha):
    """
    RDT-modified TOV equations with realistic EOS.
    
    Args:
        r: radial coordinate [cm]
        y: [P, m] = [pressure, enclosed mass]
        eos: RealisticEOS object
        alpha: RDT density-scaling parameter
    
    Returns:
        [dP/dr, dm/dr]
    """
    P, m = y
    
    if P <= 0 or r <= 0:  # Changed from 1e20 - was stopping too early!
        return [0.0, 0.0]
    
    # Get thermodynamic quantities from EOS
    epsilon = eos.get_epsilon(P)
    rho = eos.get_rho(P)
    
    if rho <= 0 or epsilon <= 0:
        return [0.0, 0.0]
    
    # Calculate RDT geometric factor
    F = F_geom(rho, alpha)
    
    # TOV equation with RDT factor
    numerator = G * (epsilon + P) * (m + 4.0*np.pi*r**3 * P/c**2)
    denominator = c**2 * r**2 * (1.0 - 2.0*G*m/(r*c**2))
    
    if denominator <= 0:
        return [0.0, 0.0]
    
    dP_dr = -F * numerator / denominator  # RDT modification
    dm_dr = 4.0*np.pi*r**2 * epsilon / c**2
    
    if not (np.isfinite(dP_dr) and np.isfinite(dm_dr)):
        return [0.0, 0.0]
    
    return [dP_dr, dm_dr]

# ============================================================================
# INTEGRATION
# ============================================================================

def integrate_NS(P_c, eos, alpha=None, r_max=30*km, debug=False):
    """
    Integrate neutron star structure.
    
    Args:
        P_c: central pressure [dyne/cm^2]
        eos: RealisticEOS object
        alpha: RDT parameter (None for standard TOV)
        r_max: maximum radius for integration [cm]
        debug: print diagnostic info
    
    Returns:
        M: total gravitational mass [Msun]
        R: stellar radius [km]
        sol: integration solution object
    """
    # Initial conditions at small radius
    r0 = 1.0e4  # 0.1 km
    rho_c = eos.get_rho(P_c)
    m0 = (4.0/3.0) * np.pi * r0**3 * rho_c
    y0 = [P_c, m0]
    
    if debug:
        print(f"  Debug: P_c = {P_c:.3e}, rho_c = {rho_c:.3e}")
        print(f"  Initial: r0 = {r0} cm, P0 = {P_c:.3e}, m0 = {m0:.3e} g")
        
        # Test derivatives at r0
        if alpha is None:
            derivs = tov_rhs_standard(r0, y0, eos)
        else:
            derivs = tov_rhs_rdt(r0, y0, eos, alpha)
        print(f"  Test derivatives at r0:")
        print(f"    dP/dr = {derivs[0]:.3e}")
        print(f"    dm/dr = {derivs[1]:.3e}")
        
        eps = eos.get_epsilon(P_c)
        metric = 1.0 - 2.0*G*m0/(r0*c**2)
        print(f"  Energy density: {eps:.3e} erg/cm^3")
        print(f"  Metric factor (1 - 2GM/(rc²)): {metric:.6f}")
        print(f"  Compactness GM/(rc^2): {G*m0/(r0*c**2):.3e}")
    
    # Event: stop when pressure reaches near-zero (surface)
    def event_surface(r, y):
        return y[0] - 1e25  # Small but non-zero threshold
    
    event_surface.terminal = True
    event_surface.direction = -1.0
    
    # Choose RHS function
    if alpha is None:
        rhs = lambda r, y: tov_rhs_standard(r, y, eos)
    else:
        rhs = lambda r, y: tov_rhs_rdt(r, y, eos, alpha)
    
    # Integrate
    sol = solve_ivp(
        rhs,
        (r0, r_max),
        y0,
        events=event_surface,
        max_step=1.0*km,
        rtol=1e-6,
        atol=1e-8,
        method='DOP853',
        dense_output=False
    )
    
    if debug:
        print(f"  Integration status: {sol.status}")
        print(f"  Integration message: {sol.message}")
        print(f"  Number of points: {len(sol.t)}")
        print(f"  Final radius: {sol.t[-1]/km:.2f} km")
        print(f"  Final P: {sol.y[0,-1]:.3e}")
        print(f"  Event triggered: {len(sol.t_events[0]) > 0}")
    
    # Check if surface was found
    if len(sol.t_events[0]) > 0:
        R = sol.t_events[0][0]
        M = np.interp(R, sol.t, sol.y[1])
        if debug:
            print(f"  Surface found: R = {R/km:.2f} km, M = {M/Msun:.4f} Msun")
    else:
        if debug:
            print(f"  Surface not found (integration reached r_max or failed)")
        return None, None, sol
    
    # Convert to standard units
    M_Msun = M / Msun
    R_km = R / km
    
    return M_Msun, R_km, sol

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    print("="*70)
    print("PAPER III PHASE 2: REALISTIC EOS TOV WITH RDT")
    print("="*70)
    
    # Load EOSs
    print("\nLoading equations of state...")
    eos_sly4 = RealisticEOS(
        '/mnt/user-data/outputs/eos_sly4.dat',
        name="SLy4"
    )
    print()
    eos_apr4 = RealisticEOS(
        '/mnt/user-data/outputs/eos_apr4.dat',
        name="APR4"
    )
    
    print("\n" + "="*70)
    print("DIAGNOSTIC TEST: Single integration with debug output")
    print("="*70)
    
    # Test with SLy4, P_c that should give ~1.4 Msun star
    P_test = 2e35  # dyne/cm²
    print(f"\nTesting SLy4 with P_c = {P_test:.2e} dyne/cm²")
    M_test, R_test, sol_test = integrate_NS(P_test, eos_sly4, alpha=None, debug=True)
    
    if M_test is not None:
        print(f"\n✓ SUCCESS: M = {M_test:.4f} Msun, R = {R_test:.2f} km")
    else:
        print(f"\n✗ FAILED: Integration did not find surface")
        print("Exiting...")
        sys.exit(1)
    
    print("\n" + "="*70)
    
    # Define scan parameters
    n_pressures = 40
    P_c_min = 5e34
    P_c_max = 5e36
    P_c_array = np.logspace(np.log10(P_c_min), np.log10(P_c_max), n_pressures)
    
    alpha_values = ['Standard', 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    
    print(f"\nScanning {n_pressures} central pressures")
    print(f"  P_c range: {P_c_min:.2e} to {P_c_max:.2e} dyne/cm^2")
    print(f"Alpha values: {alpha_values}\n")
    
    # Storage for results
    results = {eos_name: {} for eos_name in ['SLy4', 'APR4']}
    
    # Scan both EOSs
    for eos_name, eos in [('SLy4', eos_sly4), ('APR4', eos_apr4)]:
        
        print(f"\n{'='*70}")
        print(f"SCANNING: {eos_name}")
        print(f"{'='*70}\n")
        
        for alpha_label in alpha_values:
            
            alpha_val = None if alpha_label == 'Standard' else alpha_label
            
            print(f"Processing alpha = {alpha_label}")
            
            masses = []
            radii = []
            rho_centers = []
            
            success_count = 0
            
            for i, P_c in enumerate(P_c_array):
                M, R, sol = integrate_NS(P_c, eos, alpha=alpha_val)
                
                if M is not None and R is not None:
                    masses.append(M)
                    radii.append(R)
                    rho_c = eos.get_rho(P_c)
                    rho_centers.append(rho_c)
                    success_count += 1
                else:
                    masses.append(np.nan)
                    radii.append(np.nan)
                    rho_centers.append(np.nan)
                
                if (i+1) % 10 == 0:
                    print(f"  Progress: {i+1}/{n_pressures} models, {success_count} successful")
            
            print(f"  Total: {success_count}/{n_pressures} integrations succeeded\n")
            
            results[eos_name][alpha_label] = {
                'M': np.array(masses),
                'R': np.array(radii),
                'rho_c': np.array(rho_centers),
                'P_c': P_c_array
            }
    
    # Find maximum masses
    print("\n" + "="*70)
    print("TABLE: Maximum Mass Results (Phase 2 - Realistic EOS)")
    print("="*70)
    
    max_mass_table = []
    
    for eos_name in ['SLy4', 'APR4']:
        print(f"\n{eos_name}:")
        print(f"{'Alpha':<10} {'M_max (Msun)':<15} {'R (km)':<12} {'rho_c (g/cm3)':<15} {'d_eff':<10}")
        print("-"*70)
        
        for alpha_label in alpha_values:
            data = results[eos_name][alpha_label]
            masses = data['M']
            radii = data['R']
            rho_centers = data['rho_c']
            
            # Find maximum mass
            valid = ~np.isnan(masses)
            if np.any(valid):
                i_max = np.nanargmax(masses)
                M_max = masses[i_max]
                R_max = radii[i_max]
                rho_max = rho_centers[i_max]
                
                # Calculate d_eff at maximum mass
                if alpha_label == 'Standard':
                    d_eff = 3.0
                else:
                    d_eff = d_eff_spatial(rho_max, alpha_label)
                
                print(f"{str(alpha_label):<10} {M_max:<15.4f} {R_max:<12.3f} {rho_max:<15.3e} {d_eff:<10.4f}")
                
                max_mass_table.append({
                    'EOS': eos_name,
                    'alpha': alpha_label,
                    'M_max': M_max,
                    'R': R_max,
                    'rho_c': rho_max,
                    'd_eff': d_eff
                })
    
    # Day 3 Decision Gate
    print("\n" + "="*70)
    print("DAY 3 DECISION GATE: M_max > 2.1 Msun CHECK")
    print("="*70)
    
    all_pass = True
    failures = []
    
    for row in max_mass_table:
        status = "PASS" if row['M_max'] > 2.1 else "FAIL"
        if status == "FAIL":
            all_pass = False
            failures.append(f"  {row['EOS']} alpha={row['alpha']}: M_max = {row['M_max']:.4f} Msun < 2.1")
        
        label = f"{row['EOS']:<10} alpha={str(row['alpha']):<8}"
        print(f"{label} M_max = {row['M_max']:.4f} Msun ... {status}")
    
    if all_pass:
        print("\n✓ ALL MODELS PASS M_max > 2.1 Msun CONSTRAINT")
        print("✓ PROCEED TO PHASE 3 (MANUSCRIPT PREPARATION)")
    else:
        print(f"\n✗ CONSTRAINT VIOLATED FOR:")
        for failure in failures:
            print(failure)
        print("\n⚠ INVESTIGATE BEFORE PROCEEDING")
    
    print("="*70)
    
    # Save results
    print("\nSaving results...")
    
    for eos_name in ['SLy4', 'APR4']:
        filename = f'/mnt/user-data/outputs/results_phase2_{eos_name.lower()}.csv'
        
        # Combine all alpha results
        df_list = []
        for alpha_label in alpha_values:
            data = results[eos_name][alpha_label]
            df = pd.DataFrame({
                'P_c': data['P_c'],
                'M': data['M'],
                'R': data['R'],
                'rho_c': data['rho_c'],
                'alpha': str(alpha_label)
            })
            df_list.append(df)
        
        df_combined = pd.concat(df_list, ignore_index=True)
        df_combined.to_csv(filename, index=False)
        print(f"  Saved: {filename}")
    
    print("\n" + "="*70)
    print("PHASE 2 INTEGRATION COMPLETE")
    print("="*70)
