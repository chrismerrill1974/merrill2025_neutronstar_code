"""
Realistic EOS Class for Neutron Star TOV Integration

Reads tabulated EOS data and provides interpolation functions
Compatible with Phase 1 TOV code structure
"""

import numpy as np
from scipy.interpolate import interp1d
import os

class RealisticEOS:
    """
    Tabulated equation of state with interpolation
    
    Attributes:
        name: EOS name (e.g., "SLy4", "APR4")
        rho: mass density table [g/cm³]
        P: pressure table [dyne/cm²]
        epsilon: energy density table [erg/cm³]
        n_b: baryon number density table [1/cm³]
    """
    
    def __init__(self, filename, name="EOS"):
        """
        Initialize EOS from table file
        
        Args:
            filename: path to EOS table (ASCII format)
            name: descriptive name for this EOS
        """
        self.name = name
        self.filename = filename
        
        # Load table
        self._load_table()
        
        # Create interpolators
        self._create_interpolators()
        
        print(f"Loaded EOS: {self.name}")
        print(f"  Points: {len(self.rho)}")
        print(f"  Density range: {self.rho.min():.2e} to {self.rho.max():.2e} g/cm³")
        print(f"  Pressure range: {self.P.min():.2e} to {self.P.max():.2e} dyne/cm²")
    
    def _load_table(self):
        """Load EOS table from file"""
        # Read data (skip comment lines starting with #)
        data = np.loadtxt(self.filename)
        
        self.rho = data[:, 0]      # g/cm³
        self.P = data[:, 1]        # dyne/cm²
        self.epsilon = data[:, 2]  # erg/cm³
        self.n_b = data[:, 3]      # 1/cm³
        
        # Store limits for bounds checking
        self.rho_min = self.rho.min()
        self.rho_max = self.rho.max()
        self.P_min = self.P.min()
        self.P_max = self.P.max()
    
    def _create_interpolators(self):
        """
        Create interpolation functions
        
        Use log-log interpolation for better behavior across
        many orders of magnitude
        """
        # Log quantities for interpolation
        log_rho = np.log10(self.rho)
        log_P = np.log10(self.P)
        log_epsilon = np.log10(self.epsilon)
        log_nb = np.log10(self.n_b)
        
        # P(ρ) - pressure as function of density
        self._P_of_rho = interp1d(
            log_rho, log_P,
            kind='cubic',
            bounds_error=False,
            fill_value=(log_P[0], log_P[-1])
        )
        
        # ε(P) - energy density as function of pressure
        self._epsilon_of_P = interp1d(
            log_P, log_epsilon,
            kind='cubic',
            bounds_error=False,
            fill_value=(log_epsilon[0], log_epsilon[-1])
        )
        
        # n_b(P) - baryon density as function of pressure
        self._nb_of_P = interp1d(
            log_P, log_nb,
            kind='cubic',
            bounds_error=False,
            fill_value=(log_nb[0], log_nb[-1])
        )
        
        # ρ(P) - density as function of pressure (for inverse lookup)
        self._rho_of_P = interp1d(
            log_P, log_rho,
            kind='cubic',
            bounds_error=False,
            fill_value=(log_rho[0], log_rho[-1])
        )
    
    def get_pressure(self, rho):
        """
        Get pressure from density
        
        Args:
            rho: mass density [g/cm³]
        
        Returns:
            P: pressure [dyne/cm²]
        """
        if rho <= 0:
            return 0.0
        
        log_rho = np.log10(rho)
        log_P = self._P_of_rho(log_rho)
        return 10**log_P
    
    def get_epsilon(self, P):
        """
        Get energy density from pressure
        
        Args:
            P: pressure [dyne/cm²]
        
        Returns:
            epsilon: energy density [erg/cm³]
        """
        if P <= 0:
            return 0.0
        
        log_P = np.log10(P)
        log_epsilon = self._epsilon_of_P(log_P)
        return 10**log_epsilon
    
    def get_nb(self, P):
        """
        Get baryon number density from pressure
        
        Args:
            P: pressure [dyne/cm²]
        
        Returns:
            n_b: baryon number density [1/cm³]
        """
        if P <= 0:
            return 0.0
        
        log_P = np.log10(P)
        log_nb = self._nb_of_P(log_P)
        return 10**log_nb
    
    def get_rho(self, P):
        """
        Get mass density from pressure (baryon density in g/cm³)
        
        Args:
            P: pressure [dyne/cm²]
        
        Returns:
            rho: baryon mass density [g/cm³]
        """
        if P <= 0:
            return 0.0
        
        log_P = np.log10(P)
        log_rho = self._rho_of_P(log_P)
        return 10**log_rho
    
    def get_cs2(self, P):
        """
        Get speed of sound squared: c_s² = dP/dε
        
        Args:
            P: pressure [dyne/cm²]
        
        Returns:
            cs2: c_s²/c² (dimensionless)
        """
        # Numerical derivative
        dP = P * 0.01  # 1% perturbation
        
        eps1 = self.get_epsilon(P - dP)
        eps2 = self.get_epsilon(P + dP)
        
        if eps2 - eps1 > 0:
            cs2 = (2*dP) / (eps2 - eps1)
        else:
            cs2 = 0.0
        
        return cs2
    
    def summary(self):
        """Print summary of EOS properties"""
        print(f"\n{'='*70}")
        print(f"EOS: {self.name}")
        print(f"{'='*70}")
        print(f"File: {self.filename}")
        print(f"Points: {len(self.rho)}")
        print(f"\nDensity range:")
        print(f"  ρ_min = {self.rho_min:.3e} g/cm³")
        print(f"  ρ_max = {self.rho_max:.3e} g/cm³")
        print(f"\nPressure range:")
        print(f"  P_min = {self.P_min:.3e} dyne/cm²")
        print(f"  P_max = {self.P_max:.3e} dyne/cm²")
        
        # Check causality at a few points
        print(f"\nCausality check (c_s²/c² < 1):")
        for i in [len(self.P)//4, len(self.P)//2, 3*len(self.P)//4]:
            P_test = self.P[i]
            cs2 = self.get_cs2(P_test)
            status = "✓" if cs2 < 1.0 else "⚠"
            print(f"  P = {P_test:.2e}: c_s²/c² = {cs2:.4f} {status}")


# ============================================================================
# Test the class
# ============================================================================

if __name__ == "__main__":
    print("="*70)
    print("TESTING REALISTIC EOS CLASS")
    print("="*70)
    
    # Load both EOSs
    eos_sly4 = RealisticEOS(
        '/mnt/user-data/outputs/eos_sly4.dat',
        name="SLy4"
    )
    
    print()
    
    eos_apr4 = RealisticEOS(
        '/mnt/user-data/outputs/eos_apr4.dat',
        name="APR4"
    )
    
    # Test at a few densities
    print("\n" + "="*70)
    print("INTERPOLATION TESTS")
    print("="*70)
    
    test_densities = [1e14, 2.7e14, 5e14, 1e15, 2e15]
    
    for rho_test in test_densities:
        print(f"\nAt ρ = {rho_test:.2e} g/cm³:")
        
        for eos in [eos_sly4, eos_apr4]:
            P = eos.get_pressure(rho_test)
            eps = eos.get_epsilon(P)
            nb = eos.get_nb(P)
            cs2 = eos.get_cs2(P)
            
            print(f"  {eos.name:5s}: P = {P:.3e}, ε = {eps:.3e}, c_s²/c² = {cs2:.4f}")
    
    # Test inverse lookup
    print("\n" + "="*70)
    print("INVERSE LOOKUP TEST")
    print("="*70)
    
    P_test = 1e35
    print(f"\nAt P = {P_test:.2e} dyne/cm²:")
    
    for eos in [eos_sly4, eos_apr4]:
        rho = eos.get_rho(P_test)
        eps = eos.get_epsilon(P_test)
        nb = eos.get_nb(P_test)
        
        # Verify by going back
        P_check = eos.get_pressure(rho)
        
        print(f"  {eos.name}:")
        print(f"    ρ = {rho:.3e} g/cm³")
        print(f"    ε = {eps:.3e} erg/cm³")
        print(f"    n_b = {nb:.3e} 1/cm³")
        print(f"    Round-trip: P = {P_check:.3e} (error: {abs(P_check-P_test)/P_test*100:.2e}%)")
    
    # Print summaries
    eos_sly4.summary()
    eos_apr4.summary()
    
    print("\n" + "="*70)
    print("EOS CLASS READY FOR TOV INTEGRATION")
    print("="*70)
