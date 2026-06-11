#!/usr/bin/env python3
"""
Validation script for Paper III neutron star code
Checks that EOS implementation produces correct results

Version: 1.0 (Final)
Date: December 7, 2025
"""

import sys
import os
import numpy as np

# Add parent directory to path if needed
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from realistic_eos import RealisticEOS
    from paper3_phase2_realistic import integrate_NS
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Make sure you're running from the package root directory")
    sys.exit(1)

def validate_eos():
    """Validate EOS implementation against literature"""
    
    print("="*70)
    print("VALIDATING PAPER III NEUTRON STAR CODE")
    print("="*70)
    print()
    
    # Try to load corrected EOS
    eos_file = 'data/eos_sly4_dh2001_v2.dat'
    
    if not os.path.exists(eos_file):
        print(f"❌ CRITICAL: Corrected EOS file not found: {eos_file}")
        print()
        print("This means you may have the WRONG version (V2)!")
        print("Required file: eos_sly4_dh2001_v2.dat")
        print()
        return False
    
    print(f"✓ Found corrected EOS file: {eos_file}")
    print()
    
    # Load EOS
    try:
        eos = RealisticEOS(eos_file, 'SLy4-Corrected')
        print(f"✓ Loaded EOS: {eos.name}")
        print(f"  Points: {eos.npoints}")
        print(f"  Density range: {eos.rho_min:.2e} to {eos.rho_max:.2e} g/cm³")
        print()
    except Exception as e:
        print(f"❌ Failed to load EOS: {e}")
        return False
    
    # Literature values from Douchin & Haensel (2001)
    literature = {
        1.0: 13.0,   # km
        1.2: 12.3,
        1.4: 11.7,   # Canonical
        1.6: 11.5,
        1.8: 11.2,
    }
    
    # Tolerance for "acceptable" agreement
    acceptable_error = 10.0  # percent
    excellent_error = 7.0    # percent
    
    print("Testing key neutron star configurations...")
    print("-" * 70)
    print(f"{'Mass':<8} {'R_code':<10} {'R_lit':<10} {'Error':<10} {'Status':<10}")
    print(f"{'(Msun)':<8} {'(km)':<10} {'(km)':<10} {'(%)':<10} {'':<10}")
    print("-" * 70)
    
    all_passed = True
    results = []
    
    # Test at different central pressures to find target masses
    P_test = {
        1.0: 5e34,
        1.2: 8e34,
        1.4: 1e35,
        1.6: 1.5e35,
        1.8: 3e35,
    }
    
    for M_target, R_lit in literature.items():
        # Try to find configuration near target mass
        P_c = P_test[M_target]
        
        try:
            M, R, sol = integrate_NS(P_c, eos, alpha=None)
            
            if not sol.success:
                print(f"{M_target:<8.1f} {'FAILED':<10} {R_lit:<10.2f} {'N/A':<10} {'❌ FAIL':<10}")
                all_passed = False
                continue
            
            # Calculate error
            error = (R - R_lit) / R_lit * 100
            
            # Determine status
            if abs(error) < excellent_error:
                status = "✓✓ EXCELLENT"
            elif abs(error) < acceptable_error:
                status = "✓ GOOD"
            else:
                status = "⚠ HIGH ERROR"
                all_passed = False
            
            print(f"{M:<8.2f} {R:<10.2f} {R_lit:<10.2f} {error:>+9.1f} {status:<10}")
            results.append((M, R, R_lit, error))
            
        except Exception as e:
            print(f"{M_target:<8.1f} {'ERROR':<10} {R_lit:<10.2f} {'N/A':<10} {'❌ FAIL':<10}")
            print(f"  Error: {e}")
            all_passed = False
    
    print("-" * 70)
    print()
    
    # Check canonical configuration specifically
    print("Checking CANONICAL 1.4 Msun configuration...")
    print()
    
    M_14, R_14, R_lit_14, error_14 = None, None, 11.7, None
    for M, R, R_lit, error in results:
        if abs(M - 1.4) < 0.05:
            M_14, R_14, error_14 = M, R, error
            break
    
    if R_14 is not None:
        print(f"  Mass:             {M_14:.3f} Msun")
        print(f"  Radius (code):    {R_14:.2f} km")
        print(f"  Radius (lit):     {R_lit_14:.2f} km")
        print(f"  Error:            {error_14:+.1f}%")
        print()
        
        if abs(error_14) < excellent_error:
            print("  ✓✓ EXCELLENT agreement!")
        elif abs(error_14) < acceptable_error:
            print("  ✓ GOOD agreement")
        else:
            print("  ⚠ ERROR TOO HIGH")
            all_passed = False
        print()
        
        # Check if it's the V2 bug (R ~ 16 km)
        if R_14 > 15.0:
            print("  ❌ CRITICAL: Radius > 15 km suggests V2 EOS bug!")
            print("  You are using the WRONG EOS file!")
            print()
            return False
        
        # Check if it's approximately correct (V3)
        if 10.5 <= R_14 <= 12.0:
            print("  ✓ Radius in expected range (10.5-12.0 km)")
            print("  You are using the CORRECTED EOS (V3/Final)")
        print()
    else:
        print("  ❌ Failed to find 1.4 Msun configuration")
        all_passed = False
        print()
    
    # Test RDT modification
    print("Testing RDT modification (α=0.30)...")
    print()
    
    if M_14 and R_14:
        try:
            M_rdt, R_rdt, sol = integrate_NS(P_test[1.4], eos, alpha=0.30)
            
            if sol.success:
                delta_R = (R_rdt - R_14) / R_14 * 100
                print(f"  R_GR:   {R_14:.2f} km")
                print(f"  R_RDT:  {R_rdt:.2f} km")
                print(f"  ΔR/R:   {delta_R:+.2f}%")
                print()
                
                # Expected ~2% increase
                if 1.5 <= delta_R <= 3.0:
                    print("  ✓ RDT shift in expected range (+1.5% to +3%)")
                else:
                    print(f"  ⚠ RDT shift outside expected range (got {delta_R:+.2f}%)")
                print()
        except Exception as e:
            print(f"  ❌ RDT integration failed: {e}")
            print()
    
    # Final verdict
    print("="*70)
    if all_passed and R_14 is not None and 10.5 <= R_14 <= 12.0:
        print("✅ VALIDATION PASSED")
        print()
        print("Your code is using the CORRECTED EOS (V3/Final)")
        print("Results match literature within acceptable tolerances")
        print()
        return True
    else:
        print("❌ VALIDATION FAILED")
        print()
        if R_14 and R_14 > 15.0:
            print("CRITICAL: You are using the WRONG EOS file (V2)!")
            print()
            print("Required: eos_sly4_dh2001_v2.dat")
            print("DO NOT USE: eos_sly4_hp2004.dat (has bugs)")
        else:
            print("Some tests failed - check error messages above")
        print()
        return False

if __name__ == "__main__":
    success = validate_eos()
    sys.exit(0 if success else 1)
