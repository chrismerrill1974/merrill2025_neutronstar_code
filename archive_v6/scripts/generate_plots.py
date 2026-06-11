"""
Generate Publication-Quality M-R Curves for Paper III

Creates:
1. M-R curves for both EOSs (standard vs RDT)
2. Radius shift vs mass plot
3. Effective dimension plot vs radius
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set publication quality defaults
mpl.rcParams['font.size'] = 11
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

print("="*70)
print("GENERATING PUBLICATION-QUALITY PLOTS")
print("="*70)

# Load M-R data
sly4_std = np.loadtxt('/mnt/user-data/outputs/mr_curve_sly4_standard.dat')
sly4_rdt = np.loadtxt('/mnt/user-data/outputs/mr_curve_sly4_rdt.dat')
apr_std = np.loadtxt('/mnt/user-data/outputs/mr_curve_apr_standard.dat')
apr_rdt = np.loadtxt('/mnt/user-data/outputs/mr_curve_apr_rdt.dat')

# =========================================================================
# FIGURE 1: M-R Curves (Both EOSs, Standard vs RDT)
# =========================================================================
print("\nGenerating Figure 1: M-R curves...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: SLy4
ax1.plot(sly4_std[:, 1], sly4_std[:, 0], 'b-', linewidth=2, label='Standard TOV')
ax1.plot(sly4_rdt[:, 1], sly4_rdt[:, 0], 'r--', linewidth=2, label='RDT (α=0.20)')

# Mark maximum masses
i_max_std = np.argmax(sly4_std[:, 0])
i_max_rdt = np.argmax(sly4_rdt[:, 0])
ax1.plot(sly4_std[i_max_std, 1], sly4_std[i_max_std, 0], 'bo', markersize=8, 
         label=f'M$_{{max}}$ = {sly4_std[i_max_std, 0]:.3f} M$_\\odot$')
ax1.plot(sly4_rdt[i_max_rdt, 1], sly4_rdt[i_max_rdt, 0], 'ro', markersize=8,
         label=f'M$_{{max}}$ = {sly4_rdt[i_max_rdt, 0]:.3f} M$_\\odot$')

# Mark 1.4 Msun
idx_1p4_std = np.argmin(np.abs(sly4_std[:, 0] - 1.4))
idx_1p4_rdt = np.argmin(np.abs(sly4_rdt[:, 0] - 1.4))
ax1.plot(sly4_std[idx_1p4_std, 1], 1.4, 'bs', markersize=6, alpha=0.7)
ax1.plot(sly4_rdt[idx_1p4_rdt, 1], 1.4, 'rs', markersize=6, alpha=0.7)

# Add 2.1 Msun line
ax1.axhline(2.1, color='gray', linestyle=':', linewidth=1.5, alpha=0.7, 
            label='Heavy pulsar constraint')

ax1.set_xlabel('Radius (km)', fontsize=12)
ax1.set_ylabel('Mass (M$_\\odot$)', fontsize=12)
ax1.set_title('SLy4 EOS', fontsize=13, fontweight='bold')
ax1.legend(loc='lower right', fontsize=9)
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_xlim(9, 18)
ax1.set_ylim(0.5, 2.5)

# Panel B: APR
ax2.plot(apr_std[:, 1], apr_std[:, 0], 'b-', linewidth=2, label='Standard TOV')
ax2.plot(apr_rdt[:, 1], apr_rdt[:, 0], 'r--', linewidth=2, label='RDT (α=0.20)')

# Mark maximum masses
i_max_std = np.argmax(apr_std[:, 0])
i_max_rdt = np.argmax(apr_rdt[:, 0])
ax2.plot(apr_std[i_max_std, 1], apr_std[i_max_std, 0], 'bo', markersize=8,
         label=f'M$_{{max}}$ = {apr_std[i_max_std, 0]:.3f} M$_\\odot$')
ax2.plot(apr_rdt[i_max_rdt, 1], apr_rdt[i_max_rdt, 0], 'ro', markersize=8,
         label=f'M$_{{max}}$ = {apr_rdt[i_max_rdt, 0]:.3f} M$_\\odot$')

# Mark 1.4 Msun
idx_1p4_std = np.argmin(np.abs(apr_std[:, 0] - 1.4))
idx_1p4_rdt = np.argmin(np.abs(apr_rdt[:, 0] - 1.4))
ax2.plot(apr_std[idx_1p4_std, 1], 1.4, 'bs', markersize=6, alpha=0.7)
ax2.plot(apr_rdt[idx_1p4_rdt, 1], 1.4, 'rs', markersize=6, alpha=0.7)

# Add 2.1 Msun line
ax2.axhline(2.1, color='gray', linestyle=':', linewidth=1.5, alpha=0.7,
            label='Heavy pulsar constraint')

ax2.set_xlabel('Radius (km)', fontsize=12)
ax2.set_ylabel('Mass (M$_\\odot$)', fontsize=12)
ax2.set_title('APR EOS', fontsize=13, fontweight='bold')
ax2.legend(loc='lower right', fontsize=9)
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_xlim(9, 18)
ax2.set_ylim(0.5, 2.5)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/figure1_mr_curves.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/mnt/user-data/outputs/figure1_mr_curves.png', dpi=300, bbox_inches='tight')
print("  Saved: figure1_mr_curves.pdf")
print("  Saved: figure1_mr_curves.png")
plt.close()

# =========================================================================
# FIGURE 2: Radius Shift vs Mass
# =========================================================================
print("\nGenerating Figure 2: Radius shifts...")

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Calculate radius shifts
dR_sly4 = sly4_rdt[:, 1] - sly4_std[:, 1]
dR_apr = apr_rdt[:, 1] - apr_std[:, 1]

# Plot vs mass
ax.plot(sly4_std[:, 0], dR_sly4, 'b-', linewidth=2.5, label='SLy4')
ax.plot(apr_std[:, 0], dR_apr, 'r-', linewidth=2.5, label='APR')

# Mark special points
ax.axvline(1.4, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='M = 1.4 M$_\\odot$')
ax.axhline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

ax.set_xlabel('Mass (M$_\\odot$)', fontsize=12)
ax.set_ylabel('$\\Delta$R (km) = R$_{RDT}$ - R$_{std}$', fontsize=12)
ax.set_title('Radius Shift from RDT (α=0.20)', fontsize=13, fontweight='bold')
ax.legend(loc='upper left', fontsize=11)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(0.5, 2.5)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/figure2_radius_shift.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/mnt/user-data/outputs/figure2_radius_shift.png', dpi=300, bbox_inches='tight')
print("  Saved: figure2_radius_shift.pdf")
print("  Saved: figure2_radius_shift.png")
plt.close()

# =========================================================================
# FIGURE 3: Fractional Shifts vs Mass
# =========================================================================
print("\nGenerating Figure 3: Fractional shifts...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Mass fractional shift
dM_dM_sly4 = 100 * (sly4_rdt[:, 0] - sly4_std[:, 0]) / sly4_std[:, 0]
dM_dM_apr = 100 * (apr_rdt[:, 0] - apr_std[:, 0]) / apr_std[:, 0]

ax1.plot(sly4_std[:, 1], dM_dM_sly4, 'b-', linewidth=2.5, label='SLy4')
ax1.plot(apr_std[:, 1], dM_dM_apr, 'r-', linewidth=2.5, label='APR')

ax1.set_xlabel('Radius (km)', fontsize=12)
ax1.set_ylabel('$\\Delta$M/M (%)', fontsize=12)
ax1.set_title('Mass Fractional Shift', fontsize=13, fontweight='bold')
ax1.legend(loc='best', fontsize=11)
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.axhline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

# Radius fractional shift
dR_dR_sly4 = 100 * (sly4_rdt[:, 1] - sly4_std[:, 1]) / sly4_std[:, 1]
dR_dR_apr = 100 * (apr_rdt[:, 1] - apr_std[:, 1]) / apr_std[:, 1]

ax2.plot(sly4_std[:, 1], dR_dR_sly4, 'b-', linewidth=2.5, label='SLy4')
ax2.plot(apr_std[:, 1], dR_dR_apr, 'r-', linewidth=2.5, label='APR')

ax2.set_xlabel('Radius (km)', fontsize=12)
ax2.set_ylabel('$\\Delta$R/R (%)', fontsize=12)
ax2.set_title('Radius Fractional Shift', fontsize=13, fontweight='bold')
ax2.legend(loc='best', fontsize=11)
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.axhline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/figure3_fractional_shifts.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/mnt/user-data/outputs/figure3_fractional_shifts.png', dpi=300, bbox_inches='tight')
print("  Saved: figure3_fractional_shifts.pdf")
print("  Saved: figure3_fractional_shifts.png")
plt.close()

# =========================================================================
# FIGURE 4: Effective Dimension vs Radius (for one EOS)
# =========================================================================
print("\nGenerating Figure 4: Effective dimension profile...")

# Load one of the solution profiles to get d_eff vs radius
# We'll need to recalculate this from the TOV solution
# For now, create a representative profile

# Import what we need
import sys
sys.path.insert(0, '/home/claude')
sys.path.insert(0, '/mnt/user-data/outputs')
from realistic_eos import RealisticEOS
from paper3_phase2_realistic import integrate_NS

# Generate one representative profile
eos_sly4 = RealisticEOS('/mnt/user-data/outputs/eos_sly4_hp2004.dat', 'SLy4')
P_c = 2.3e35  # Near maximum mass
M, R, sol = integrate_NS(P_c, eos_sly4, alpha=0.20, debug=False)

if sol is not None:
    # Extract profile
    r = sol.t / 1e5  # Convert to km
    P_profile = sol.y[0, :]  # Pressure
    
    # Get density from pressure using EOS
    rho_profile = np.array([eos_sly4.get_rho(P) if P > 0 else 0 for P in P_profile])
    
    # Calculate d_eff
    rho_0 = 150.0  # g/cm³
    A = 0.0334
    alpha = 0.20
    
    Omega_spatial = 1 - A * (rho_profile/rho_0)**alpha / (1 + (rho_profile/rho_0)**alpha)
    d_eff = 3 * Omega_spatial
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    ax.plot(r, d_eff, 'b-', linewidth=2.5)
    ax.axhline(3, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, 
               label='Standard 3D')
    
    ax.set_xlabel('Radius (km)', fontsize=12)
    ax.set_ylabel('Effective Spatial Dimension $d_{eff}$', fontsize=12)
    ax.set_title('Dimensional Opening Profile (SLy4, α=0.20)', fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_ylim(2.8, 3.05)
    
    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/figure4_deff_profile.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('/mnt/user-data/outputs/figure4_deff_profile.png', dpi=300, bbox_inches='tight')
    print("  Saved: figure4_deff_profile.pdf")
    print("  Saved: figure4_deff_profile.png")
    plt.close()

print("\n" + "="*70)
print("PLOT GENERATION COMPLETE")
print("="*70)
print("\nGenerated 4 figures:")
print("  1. M-R curves (both EOSs, standard vs RDT)")
print("  2. Radius shift vs mass")
print("  3. Fractional shifts vs radius")
print("  4. Effective dimension profile")
print("\nAll figures saved in PDF and PNG formats.")
