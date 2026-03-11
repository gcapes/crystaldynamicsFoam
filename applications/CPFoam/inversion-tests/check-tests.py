#!/usr/bin/env python3
"""
check-tests.py — Validates output from the 4 CPFoam inversion test cases.

must have exit(0); in createFields after inversion with verbose = true to be run effectively. 

Run:
    cd inversion-tests && ./Allrun_tests
    python3 check-tests.py
"""

import re
import sys
import os

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

# helpers

GREEN = "\033[32m"
RED   = "\033[31m"
BOLD  = "\033[1m"
RESET = "\033[0m"

passed = 0
failed = 0


def ok(msg):
    global passed
    passed += 1
    print(f"  {GREEN}PASS{RESET}  {msg}")


def fail(msg, detail=""):
    global failed
    failed += 1
    line = f"  {RED}FAIL{RESET}  {msg}"
    if detail:
        line += f"\n        → {detail}"
    print(line)


def read_log(test_name):
    path = os.path.join(TESTS_DIR, test_name, "log.CPFoam")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return f.read()


def check_symmetry_diagnostics(log, label):
    """All three residuals must be reported as 0."""
    patterns = [
        r"max \|C\(i,j,k,l\) - C\(j,i,k,l\)\| = (\S+)",
        r"max \|C\(i,j,k,l\) - C\(i,j,l,k\)\| = (\S+)",
        r"max \|C\(i,j,k,l\) - C\(k,l,i,j\)\| = (\S+)",
    ]
    names = ["minor(ij)", "minor(kl)", "major"]
    all_ok = True
    for pat, name in zip(patterns, names):
        m = re.search(pat, log)
        if not m:
            fail(f"{label}: symmetry diagnostic '{name}' not found in log")
            all_ok = False
        else:
            val = float(m.group(1))
            if val == 0.0:
                ok(f"{label}: {name} residual = 0")
            else:
                fail(f"{label}: {name} residual = {val} (expected 0 after castSymmetry)")
                all_ok = False
    return all_ok


def parse_mandel_diagonal(log):
    """Return the 6 diagonal values from 'M_inv (6x6)' block."""
    # Find the matrix block
    m = re.search(r"M_inv \(6x6\).*?:\n((?:.*\n){6})", log)
    if not m:
        return None
    rows = m.group(1).strip().split("\n")
    diag = []
    for i, row in enumerate(rows):
        vals = row.split()
        if len(vals) >= i + 1:
            diag.append(float(vals[i]))
    return diag if len(diag) == 6 else None


def check_complete(log, label):
    if "...complete" in log:
        ok(f"{label}: invertRankFour pipeline completed")
        return True
    else:
        fail(f"{label}: '...complete' not found — pipeline did not finish")
        return False


def approx_equal(a, b, rtol=1e-2):
    """Relative tolerance check."""
    if b == 0:
        return abs(a) < 1e-20
    return abs(a - b) / abs(b) < rtol


#test-1: Tom's original cubic case, C_max~1e10

def check_test1():
    print(f"\n{BOLD}inversion-test-1{RESET}  (cubic FCC, C_max~1e10, phi1=0.26)")
    log = read_log("inversion-test-1")
    if log is None:
        fail("test-1: log.CPFoam not found - did you run ./Allrun_tests?")
        return

    check_complete(log, "test-1")
    check_symmetry_diagnostics(log, "test-1")

    diag = parse_mandel_diagonal(log)
    if diag is None:
        fail("test-1: could not parse M_inv diagonal")
    else:
        # Expect S_M(0,0) ~ 1.49e-10 (cubic, rotated)
        if approx_equal(diag[0], 1.49e-10, rtol=0.05):
            ok(f"test-1: M_inv(0,0) = {diag[0]:.6e}  (expected ~1.49e-10)")
        else:
            fail(f"test-1: M_inv(0,0) = {diag[0]:.6e}  (expected ~1.49e-10)")


#test-2: cubic FCC, C_max~1e11 (10× stiffer)

def check_test2():
    print(f"\n{BOLD}inversion-test-2{RESET}  (cubic FCC, C_max~1e11 - 10x test-1, phi1=0.26)")
    log = read_log("inversion-test-2")
    if log is None:
        fail("test-2: log.CPFoam not found - did you run ./Allrun_tests?")
        return

    check_complete(log, "test-2")
    check_symmetry_diagnostics(log, "test-2")

    diag = parse_mandel_diagonal(log)
    if diag is None:
        fail("test-2: could not parse M_inv diagonal")
    else:
        # 10× stiffer → compliance 10× smaller: ~1.49e-11
        if approx_equal(diag[0], 1.49e-11, rtol=0.05):
            ok(f"test-2: M_inv(0,0) = {diag[0]:.6e}  (expected ~1.49e-11, 10x less than test-1)")
        else:
            fail(f"test-2: M_inv(0,0) = {diag[0]:.6e}  (expected ~1.49e-11)")

        log1 = read_log("inversion-test-1")
        if log1:
            diag1 = parse_mandel_diagonal(log1)
            if diag1:
                ratio = diag1[0] / diag[0]
                if approx_equal(ratio, 10.0, rtol=0.02):
                    ok(f"test-2: compliance ratio test-1/test-2 = {ratio:.4f}  (expected 10.0)")
                else:
                    fail(f"test-2: compliance ratio test-1/test-2 = {ratio:.4f}  (expected 10.0)")


#test-3: isotropic E=200GPa, nu=0.3, analytical compliance known 

def check_test3():
    print(f"\n{BOLD}inversion-test-3{RESET}  (isotropic E=200 GPa, nu=0.3, analytical compliance)")
    log = read_log("inversion-test-3")
    if log is None:
        fail("test-3: log.CPFoam not found - did you run ./Allrun_tests?")
        return

    check_complete(log, "test-3")
    check_symmetry_diagnostics(log, "test-3")

    # Analytical Mandel compliance (isotropic, unrotated):
    #   S_M(0,0) = 1/E = 5.000e-12
    #   S_M(0,1) = -nu/E = -1.500e-12
    #   S_M(3,3) = (1+nu)/E * 2  (Mandel factor) = 2*(1.3)/200e9 = 1/E * 2*(1+nu)
    #            but stored as 1/(2G) in Mandel = 1/(E/(1+nu)) = (1+nu)/E
    #            Mandel shear: 2*(1+nu)/E = 2*1.3/200e9 = 13e-12 / 2 = 6.5e-12
    # (The Mandel 6x6 shear block stores 2*(1+nu)/E after the sqrt(2) factors cancel.)
    diag = parse_mandel_diagonal(log)
    if diag is None:
        fail("test-3: could not parse M_inv diagonal")
        return

    # Normal diagonal
    for i, expected in zip([0, 1, 2], [5.0e-12, 5.0e-12, 5.0e-12]):
        if approx_equal(diag[i], expected, rtol=0.01):
            ok(f"test-3: M_inv({i},{i}) = {diag[i]:.6e}  (expected {expected:.3e} = 1/E)")
        else:
            fail(f"test-3: M_inv({i},{i}) = {diag[i]:.6e}  (expected {expected:.3e})")

    # Shear diagonal (Mandel)
    for i, expected in zip([3, 4, 5], [6.5e-12, 6.5e-12, 6.5e-12]):
        if approx_equal(diag[i], expected, rtol=0.01):
            ok(f"test-3: M_inv({i},{i}) = {diag[i]:.6e}  (expected {expected:.3e} = 2(1+nu)/E)")
        else:
            fail(f"test-3: M_inv({i},{i}) = {diag[i]:.6e}  (expected {expected:.3e})")

    # Off-diagonal: parse row 0 to get S_M(0,1)
    m = re.search(r"M_inv \(6x6\).*?:\n(.*)\n", log)
    if m:
        row0_vals = m.group(1).split()
        if len(row0_vals) >= 2:
            s01 = float(row0_vals[1])
            if approx_equal(s01, -1.5e-12, rtol=0.01):
                ok(f"test-3: M_inv(0,1) = {s01:.6e}  (expected -1.5e-12 = -nu/E)")
            else:
                fail(f"test-3: M_inv(0,1) = {s01:.6e}  (expected -1.5e-12)")


# test-4: intentionally asymmetric input - warnings + symmetrisation

def check_test4():
    print(f"\n{BOLD}inversion-test-4{RESET}  (intentionally asymmetric input - warning + castSymmetry pipeline)")
    log = read_log("inversion-test-4")
    if log is None:
        fail("test-4: log.CPFoam not found - did you run ./Allrun_tests?")
        return

    # Warnings should be present (non-symmetric pairs flagged at init)
    warning_count = log.count("--> FOAM Warning :")
    if warning_count >= 6:
        ok(f"test-4: {warning_count} FOAM Warning(s) fired for asymmetric pairs (expected >=6)")
    else:
        fail(f"test-4: only {warning_count} FOAM Warning(s) (expected ≥6 for asymmetric pairs)")

    # No FatalError
    if "FOAM FATAL ERROR" in log or "FatalError" in log:
        fail("test-4: FatalError found - castSymmetry did not fix asymmetry before fourthToMandel")
    else:
        ok("test-4: no FatalError: assertSymmetry passed after castSymmetry")

    check_complete(log, "test-4")
    check_symmetry_diagnostics(log, "test-4")

    # Symmetrised C0011 ~ avg(9.20e10, 8.50e10) = 8.85e10
    # Compliance diagonal should be slightly lower than test-2 (averaged-down stiffness)
    diag = parse_mandel_diagonal(log)
    if diag:
        ok(f"test-4: M_inv diagonal parsed: S_M(0,0) = {diag[0]:.6e}")


# main 

if __name__ == "__main__":
    print(f"{BOLD}CPFoam inversion test checker{RESET}")
    print("=" * 50)

    check_test1()
    check_test2()
    check_test3()
    check_test4()

    print()
    print("=" * 50)
    total = passed + failed
    colour = GREEN if failed == 0 else RED
    print(f"{colour}{BOLD}{passed}/{total} checks passed{RESET}", end="")
    if failed:
        print(f"  ({failed} failed)")
    else:
        print()

    sys.exit(0 if failed == 0 else 1)
