#!/usr/bin/env python3
"""
plotWaveSpeed.py -- crystaldynamicsFoam elastic wave-speed validation postprocessor.

Reads the OpenFOAM `probes` function-object output of the shear- and
pressure-wave runs, detects the first arrival of the wave front at each
probe, fits arrival time against probe position to obtain the numerical
propagation speed, and compares it with the analytical value derived from
constant/GrainProperties:

    shear        c_s = sqrt(C44 / rho)   [C44 read from keyword C0101]
    longitudinal c_p = sqrt(C11 / rho)   [C11 read from keyword C0000]

(uniaxial-strain plane waves along a <100> direction of the crystal frame,
which for the distributed isotropic material holds in every direction).

Usage:
    plotWaveSpeed.py SHEAR_CASE PRESSURE_CASE [-o waveSpeed.png] [--tol 0.02]

Exit status 0 if both measured speeds are within `tol` (relative) of the
analytical values, 1 otherwise -- so the script doubles as a regression
test inside Allrun.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------------
# GrainProperties parsing
# ----------------------------------------------------------------------------
def read_grain_properties(case):
    """Return (C11, C44, rho0) from constant/GrainProperties."""
    path = os.path.join(case, "constant", "GrainProperties")
    with open(path) as fh:
        text = fh.read()
    # strip // and /* */ comments
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    text = re.sub(r"//[^\n]*", " ", text)

    def grab(key):
        m = re.search(r"\b%s\s+([0-9.eE+\-]+)\s*;" % key, text)
        if not m:
            raise RuntimeError("keyword '%s' not found in %s" % (key, path))
        return float(m.group(1))

    C11, C12, C44, rho0 = (grab(k) for k in ("C0000", "C0011", "C0101",
                                             "rho_0"))
    # isotropy sanity check: C11 - C12 = 2 C44
    dev = abs(C11 - C12 - 2.0 * C44) / C11
    if dev > 1e-6:
        print("WARNING: stiffness is not isotropic "
              "(|C11-C12-2*C44|/C11 = %.3g); the analytic speeds below are "
              "only valid for propagation along a cube axis of the crystal "
              "frame." % dev)
    return C11, C44, rho0


# ----------------------------------------------------------------------------
# Probe file parsing
# ----------------------------------------------------------------------------
def read_probes(case, field="U"):
    """Return (positions[nP,3], t[nT], data[nT,nP,3]) from the probes FO."""
    pattern = os.path.join(case, "postProcessing", "probes", "*", field)
    files = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError(
            "no probe data found under %s\n"
            "(did the run finish? is the 'probes' function object in "
            "system/controlDict?)" % pattern)
    path = files[0]

    positions, times, rows = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                m = re.match(r"#\s*Probe\s+\d+\s+\(([^)]*)\)", line)
                if m:
                    positions.append([float(v) for v in m.group(1).split()])
                continue
            vecs = re.findall(r"\(([^)]*)\)", line)
            if not vecs:
                continue
            t = float(line.split("(")[0])
            row = [[float(v) for v in vec.split()] for vec in vecs]
            times.append(t)
            rows.append(row)

    positions = np.array(positions)                # (nP, 3)
    t = np.array(times)                            # (nT,)
    data = np.array(rows)                          # (nT, nP, 3)
    if data.shape[1] != positions.shape[0]:
        raise RuntimeError("probe count mismatch in %s" % path)
    return positions, t, data


# ----------------------------------------------------------------------------
# Wave-speed measurement
# ----------------------------------------------------------------------------
def measure_speed(positions, t, data, frac=0.10):
    """First-arrival regression.

    Returns (speed, axis, comp, s[nT,nP], arrivals, used_pos, fit_coeffs).
    """
    # propagation axis: the coordinate with the largest spread across probes
    axis = int(np.argmax(positions.max(axis=0) - positions.min(axis=0)))
    x = positions[:, axis]

    # signal component: the one carrying the wave (largest peak magnitude)
    comp = int(np.argmax(np.abs(data).max(axis=(0, 1))))
    s = data[:, :, comp]                           # (nT, nP)

    thr = frac * np.abs(s).max()
    arrivals = np.full(x.shape, np.nan)
    for p in range(x.size):
        idx = np.nonzero(np.abs(s[:, p]) >= thr)[0]
        if idx.size == 0:
            continue
        i = idx[0]
        if i == 0:
            arrivals[p] = t[0]
        else:
            # linear interpolation of the threshold crossing
            f0, f1 = np.abs(s[i - 1, p]), np.abs(s[i, p])
            w = (thr - f0) / max(f1 - f0, 1e-300)
            arrivals[p] = t[i - 1] + w * (t[i] - t[i - 1])

    ok = np.isfinite(arrivals)
    if ok.sum() < 3:
        raise RuntimeError("wave front reached fewer than 3 probes; "
                           "increase endTime or check the run")

    coeffs = np.polyfit(x[ok], arrivals[ok], 1)    # t = a*x + b
    speed = 1.0 / abs(coeffs[0])
    return speed, axis, comp, s, arrivals, ok, coeffs


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def plot_case(axXT, axWF, name, positions, t, s, arrivals, ok, coeffs,
              axis, comp, c_num, c_ana):
    x = positions[:, axis]
    xu = x * 1e6                                   # um
    tu = 1e9                                       # s -> ns

    # --- x-t diagram ---
    axXT.plot(xu[ok], arrivals[ok] * tu, "o", ms=5, label="probe arrivals")
    xf = np.linspace(x.min(), x.max(), 2)
    axXT.plot(xf * 1e6, np.polyval(coeffs, xf) * tu, "-",
              label="fit: %.1f m/s" % c_num)
    # analytic line anchored on the fitted arrival at the first-hit probe
    x0 = x[ok][np.argmin(arrivals[ok])]
    t0 = np.polyval(coeffs, x0)
    sgn = np.sign(coeffs[0])
    axXT.plot(xf * 1e6, (t0 + sgn * (xf - x0) / c_ana) * tu, "--k",
              label="analytic: %.1f m/s" % c_ana)
    err = 100.0 * (c_num - c_ana) / c_ana
    axXT.set_title("%s wave   error = %+.2f%%" % (name, err))
    axXT.set_xlabel(r"probe position [$\mu$m]")
    axXT.set_ylabel("first arrival [ns]")
    axXT.legend(fontsize=8)
    axXT.grid(alpha=0.3)

    # --- waterfall of the probe signals ---
    amp = np.abs(s).max()
    for p in range(x.size):
        off = p * 1.5
        axWF.plot(t * tu, s[:, p] / amp + off, lw=0.8)
        if ok[p]:
            axWF.plot(arrivals[p] * tu,
                      np.interp(arrivals[p], t, s[:, p]) / amp + off,
                      "kx", ms=6)
    axWF.set_xlabel("time [ns]")
    axWF.set_ylabel("U$_%s$ / max (offset per probe)" % "xyz"[comp])
    axWF.set_title("%s wave: probe signals (x = arrival)" % name)
    axWF.grid(alpha=0.3)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("shear_case", help="case directory of the shear run")
    ap.add_argument("pressure_case", help="case directory of the pressure run")
    ap.add_argument("-o", "--output", default="waveSpeedValidation.png")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="relative pass tolerance (default 0.02)")
    ap.add_argument("--frac", type=float, default=0.10,
                    help="arrival threshold as fraction of peak (default 0.1)")
    args = ap.parse_args()

    fig, axs = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

    summary, ok_all = [], True
    for col, (name, case) in enumerate(
            [("shear", args.shear_case), ("pressure", args.pressure_case)]):
        C11, C44, rho0 = read_grain_properties(case)
        c_ana = np.sqrt((C44 if name == "shear" else C11) / rho0)

        positions, t, data = read_probes(case)
        c_num, axis, comp, s, arr, okp, coeffs = measure_speed(
            positions, t, data, args.frac)

        # polarisation cross-check: the shear run must carry a TRANSVERSE
        # signal (component != propagation axis), the pressure run a
        # LONGITUDINAL one.  A mismatch means the intended wave was never
        # launched (e.g. a pulse polarised along an 'empty' mesh direction,
        # whose velocity component OpenFOAM excludes from the solve) and the
        # fit has locked onto the wrong wave -- fail loudly.
        want_longitudinal = (name == "pressure")
        is_longitudinal = (comp == axis)
        if is_longitudinal != want_longitudinal:
            print("ERROR [%s]: dominant signal is U_%s but the propagation "
                  "axis is %s -> the run carries a %s wave, not a %s wave.\n"
                  "       Check the pulse polarisation in 0/U (it must not "
                  "point along an 'empty' mesh direction) and the tangential "
                  "boundary conditions."
                  % (name, "xyz"[comp], "xyz"[axis],
                     "longitudinal" if is_longitudinal else "transverse",
                     "longitudinal" if want_longitudinal else "transverse"))
            ok_all = False

        err = (c_num - c_ana) / c_ana
        ok_all &= abs(err) <= args.tol
        summary.append((name, c_num, c_ana, err))

        plot_case(axs[0, col], axs[1, col], name, positions, t, s, arr,
                  okp, coeffs, axis, comp, c_num, c_ana)

    fig.suptitle("crystaldynamicsFoam elastic wave-speed validation", fontsize=13)
    fig.savefig(args.output, dpi=160)

    print("")
    print("%-10s %14s %14s %10s" % ("wave", "numerical", "analytic", "error"))
    for name, c_num, c_ana, err in summary:
        print("%-10s %11.1f m/s %11.1f m/s %+9.2f%%"
              % (name, c_num, c_ana, 100.0 * err))
    print("\nfigure written to %s" % args.output)

    if ok_all:
        print("PASS: both speeds within %.1f%%" % (100.0 * args.tol))
        return 0
    print("FAIL: at least one speed outside the %.1f%% tolerance"
          % (100.0 * args.tol))
    return 1


if __name__ == "__main__":
    sys.exit(main())
