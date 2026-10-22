#!/usr/bin/env python3
# Planar-impact Hugoniot validation: fits Us(up) of every run under runs/<mode>/up<value>/ and compares the Mie-Gruneisen mode with the input (c0, s) and the linear-elastic mode with p = K*mu.

import os
import re
import math
import sys

# case constants (must match caseTemplate)
L_BAR   = 16.0e-3        # bar length [m]
N_CELLS = 800           # cells along x
DX      = L_BAR / N_CELLS
XC      = [(i + 0.5) * DX for i in range(N_CELLS)]   # uniform-mesh centres

RHO0    = 8930.0        # [kg/m3]
C0_IN   = 3940.0        # dictionary c0        [m/s]
S_IN    = 1.489         # dictionary s_UsUp    [-]
C11, C12 = 168.4e9, 121.4e9
K_BULK  = (C11 + 2.0 * C12) / 3.0          # 137.07 GPa
CK      = math.sqrt(K_BULK / RHO0)          # linear-elastic closure sound speed

FIT_T_MIN_FRAC = 0.35   # discard the formation transient before this fraction of the latest time
PLATEAU_LO, PLATEAU_HI = 0.15, 0.50         # plateau window (fraction of x_f)

TOL_S  = 0.12           # PASS tolerance on fitted s   (MieGruneisen)
TOL_C0 = 80.0           # PASS tolerance on fitted c0  [m/s]

RUNS_DIR = "runs"

# OpenFOAM ASCII field parsing

def _matched_block(text, start):
    """Return substring inside the parenthesis block opening at text[start]=='('."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
    raise ValueError("unbalanced parentheses in field file")


def read_internal_field(fname, ncomp):
    """Parse the internalField of an ASCII volScalar/volVectorField into floats or [x,y,z] lists."""
    with open(fname) as fh:
        txt = fh.read()

    key = txt.find("internalField")
    if key < 0:
        raise ValueError(f"{fname}: no internalField")
    sub = txt[key:]

    m = re.match(r"internalField\s+(uniform|nonuniform)", sub)
    if not m:
        raise ValueError(f"{fname}: unparsable internalField")

    if m.group(1) == "uniform":
        if ncomp == 1:
            v = re.match(
                r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", sub)
            return [float(v.group(1))] * N_CELLS
        v = re.match(r"internalField\s+uniform\s+\(([^)]*)\)\s*;", sub)
        comps = [float(t) for t in v.group(1).split()]
        return [comps[:] for _ in range(N_CELLS)]

    # nonuniform: "... List<type>\n N \n ( ... )"
    mn = re.match(
        r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*", sub)
    if not mn:
        raise ValueError(f"{fname}: unparsable nonuniform header")
    n = int(mn.group(1))
    popen = sub.find("(", mn.end())
    body = _matched_block(sub, popen)

    if ncomp == 1:
        vals = [float(t) for t in body.split()]
        if len(vals) != n:
            raise ValueError(f"{fname}: expected {n} scalars, got {len(vals)}")
        return vals

    tuples = re.findall(r"\(([^)]+)\)", body)
    if len(tuples) != n:
        raise ValueError(f"{fname}: expected {n} vectors, got {len(tuples)}")
    return [[float(t) for t in tup.split()] for tup in tuples]


# analysis helpers

def front_position(ux, threshold):
    """Frontmost half-amplitude crossing of the u_x profile, or None."""
    idx = [i for i, v in enumerate(ux) if v >= threshold]
    if not idx:
        return None
    i0 = max(idx)
    if i0 >= N_CELLS - 1:
        return None                      # front has reached the far end
    u0, u1 = ux[i0], ux[i0 + 1]
    if u0 == u1:
        return XC[i0]
    frac = (u0 - threshold) / (u0 - u1)
    return XC[i0] + frac * DX


def lstsq_line(xs, ys):
    """Least-squares y = a + b x. Returns (a, b)."""
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    return a, b


def us_legacy(up):
    """Shock speed of the linear-elastic closure p = K*mu: Us = cK sqrt(1+mu), up = cK mu/sqrt(1+mu) (bisection)."""
    f = lambda mu: CK * mu / math.sqrt(1.0 + mu) - up
    lo, hi = 0.0, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    mu = 0.5 * (lo + hi)
    return CK * math.sqrt(1.0 + mu)


EXP_DATA_FILE = "cu_hugoniot_experimental.dat"

def read_experimental(path):
    """Parse (up, Us, source) rows; values < 20 are km/s and are converted to m/s."""
    pts = []
    if not os.path.isfile(path):
        return pts
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                up, Us = float(parts[0]), float(parts[1])
            except (ValueError, IndexError):
                continue
            src = " ".join(parts[2:]) if len(parts) > 2 else "experiment"
            if Us < 20.0:
                up, Us = up * 1000.0, Us * 1000.0
            pts.append((up, Us, src))
    return pts


def time_dirs(case):
    out = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if t > 0.0:
            out.append((t, os.path.join(case, d)))
    return sorted(out)


def analyse_case(case, up_nominal):
    """Return dict with Us fit, plateau up, plateau p -- or None."""
    thr = 0.5 * up_nominal
    traj = []
    for t, d in time_dirs(case):
        ufile = os.path.join(d, "U")
        if not os.path.isfile(ufile):
            continue
        try:
            ux = [v[0] for v in read_internal_field(ufile, 3)]
        except Exception as err:
            print(f"    warning: {ufile}: {err}")
            continue
        xf = front_position(ux, thr)
        if xf is not None:
            traj.append((t, xf, ux, d))
    if len(traj) < 4:
        return None

    t_last = traj[-1][0]
    fit_pts = [(t, xf) for (t, xf, _, _) in traj if t >= FIT_T_MIN_FRAC * t_last]
    if len(fit_pts) < 3:
        return None
    _, Us = lstsq_line([p[0] for p in fit_pts], [p[1] for p in fit_pts])

    # post-shock plateau at the last time: particle velocity and pressure
    t, xf, ux, d = traj[-1]
    win = [i for i, x in enumerate(XC) if PLATEAU_LO * xf <= x <= PLATEAU_HI * xf]
    up_meas = sum(ux[i] for i in win) / len(win)

    p_meas = None
    pfile = os.path.join(d, "p")
    if os.path.isfile(pfile):
        try:
            pf = read_internal_field(pfile, 1)
            p_meas = sum(pf[i] for i in win) / len(win)
        except Exception:
            pass

    return dict(Us=Us, up=up_meas, p=p_meas, traj=traj, up_nominal=up_nominal)


# main

def main():
    if not os.path.isdir(RUNS_DIR):
        sys.exit("No runs/ directory found -- run ./Allrun first.")

    modes = sorted(
        d for d in os.listdir(RUNS_DIR)
        if os.path.isdir(os.path.join(RUNS_DIR, d)))
    if not modes:
        sys.exit("runs/ is empty -- run ./Allrun first.")

    results = {}
    for mode in modes:
        mdir = os.path.join(RUNS_DIR, mode)
        cases = sorted(
            (c for c in os.listdir(mdir) if c.startswith("up")),
            key=lambda c: float(c[2:]))
        rows = []
        print(f"\n=== mode: {mode} ===")
        for c in cases:
            upn = float(c[2:])
            res = analyse_case(os.path.join(mdir, c), upn)
            if res is None:
                print(f"  up={upn:6.0f}  -- insufficient data (run failed "
                      "or front not formed); skipped")
                continue
            pj = RHO0 * res["Us"] * res["up"]
            pline = (f"  p={res['p']/1e9:7.3f} GPa "
                     f"(rho0*Us*up = {pj/1e9:7.3f} GPa)"
                     if res["p"] is not None else "")
            print(f"  up={res['up']:7.1f} m/s   Us={res['Us']:7.1f} m/s{pline}")
            rows.append(res)
        if len(rows) >= 2:
            c0f, sf = lstsq_line([r["up"] for r in rows],
                                 [r["Us"] for r in rows])
            results[mode] = (c0f, sf, rows)

    exp_pts = read_experimental(
        sys.argv[1] if len(sys.argv) > 1 else EXP_DATA_FILE)
    if exp_pts:
        srcs = {}
        for up, Us, s0 in exp_pts:
            srcs[s0] = srcs.get(s0, 0) + 1
        print(f"\nExperimental overlay: {len(exp_pts)} points from "
              + ", ".join(f"{k} ({n})" for k, n in sorted(srcs.items())))
        r_mg = [Us - (C0_IN + S_IN * up) for up, Us, _ in exp_pts]
        r_lg = [Us - us_legacy(up) for up, Us, _ in exp_pts]
        rms = lambda r: math.sqrt(sum(x * x for x in r) / len(r))
        mean = lambda r: sum(r) / len(r)
        print(f"  residuals vs MG input line : rms = {rms(r_mg):6.1f} m/s, "
              f"bias = {mean(r_mg):+6.1f} m/s")
        print(f"  residuals vs linear-elastic p=K*mu : rms = {rms(r_lg):6.1f} m/s, "
              f"bias = {mean(r_lg):+6.1f} m/s")
        print("  (the data confirming MG and rejecting the linear-elastic closure "
              "is the validation claim)")
    else:
        print("\n(no experimental overlay: populate "
              "cu_hugoniot_experimental.dat -- sources listed in its "
              "header -- to add measured shots to the figure)")

    print("\n" + "=" * 66)
    print("Fitted Hugoniots:  Us = c0_fit + s_fit * up")
    print("=" * 66)
    print(f"  dictionary input      : c0 = {C0_IN:7.1f} m/s   s = {S_IN:6.3f}")
    print(f"  linear-elastic closure p=K*mu : c0 = {CK:7.1f} m/s   s ~ 0.5 "
          "(exact curve, not a line)")
    print("-" * 66)

    passed = None
    for mode, (c0f, sf, rows) in results.items():
        print(f"  {mode:14s}        : c0 = {c0f:7.1f} m/s   s = {sf:6.3f}")
        if mode == "MieGruneisen":
            ok = abs(sf - S_IN) < TOL_S and abs(c0f - C0_IN) < TOL_C0
            passed = ok
            print(f"    -> |s - {S_IN}| = {abs(sf - S_IN):.3f} "
                  f"(tol {TOL_S}),  |c0 - {C0_IN:.0f}| = "
                  f"{abs(c0f - C0_IN):.1f} m/s (tol {TOL_C0:.0f})")
            print(f"    -> {'PASS' if ok else 'FAIL'}: EOS mode "
                  f"{'recovers' if ok else 'does NOT recover'} the input "
                  "Hugoniot")
        elif mode == "linearElastic":
            exp = [us_legacy(r["up"]) for r in rows]
            rms = math.sqrt(sum((r["Us"] - e) ** 2
                            for r, e in zip(rows, exp)) / len(rows))
            print(f"    -> rms deviation from derived p=K*mu Hugoniot: "
                  f"{rms:6.1f} m/s")
            print("    -> effective slope ~0.5 is the DEFECT the EOS "
                  "extension repairs (real Cu: s = 1.489)")

    # optional figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        up_hi = 1.05 * max(
            [900.0]
            + [r["up"] for _, (_, _, rows) in results.items() for r in rows]
            + [up for up, _, _ in exp_pts])
        ups = [i * up_hi / 180.0 for i in range(181)]
        mg_line = [C0_IN + S_IN * u for u in ups]
        # +/-0.9% band: fractional scatter of the gas-gun data about linear fits (Mitchell & Nellis 1981)
        ax.fill_between(ups, [v * 0.991 for v in mg_line],
                        [v * 1.009 for v in mg_line],
                        color="0.85", zorder=0,
                        label="$\\pm$0.9% band (gas-gun data scatter,\n"
                              "Mitchell & Nellis 1981)")
        #ax.plot(ups, mg_line, "k-",
                #label=f"Mie-Gruneisen input: $U_s = {C0_IN:.0f} + "
                     # f"{S_IN}\\,u_p$")
        #ax.plot(ups, [us_legacy(u) for u in ups], "k--",
                #label="legacy closure $p=K\\mu$ (eff. $s\\approx0.5$)")
        if exp_pts:
            mkcycle = ["^", "v", "D", "P", "X", "*", "<", ">"]
            for k, s0 in enumerate(sorted({p[2] for p in exp_pts})):
                sel = [(up, Us) for up, Us, s in exp_pts if s == s0]
                ax.plot([u for u, _ in sel], [v for _, v in sel],
                        mkcycle[k % len(mkcycle)], ms=4.5, ls="none",
                        color="0.35", mfc="none", mew=1.0, label=s0)
        marks = {"MieGruneisen": ("o", "tab:red"),
                 "linearElastic": ("s", "tab:blue")}
        for mode, (c0f, sf, rows) in results.items():
            mk, col = marks.get(mode, ("^", "tab:green"))
            ax.plot([r["up"] for r in rows], [r["Us"] for r in rows],
                    mk, color=col, ms=7, mfc="none", mew=1.6,
                    label=f"{mode}: fit $c_0={c0f:.0f}$, $s={sf:.2f}$")
        ax.set_xlabel("particle velocity  $u_p$  [m/s]")
        ax.set_ylabel("shock velocity  $U_s$  [m/s]")
        ax.set_title("crystaldynamicsFoam planar-impact Hugoniot: copper")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig("hugoniot_UsUp.png", dpi=180)
        print("\nwrote hugoniot_UsUp.png")

        for mode, (_, _, rows) in results.items():
            fig2, ax2 = plt.subplots(figsize=(6.4, 4.2))
            r = rows[-1]
            for (t, xf, ux, _) in r["traj"][::4]:
                ax2.plot([x * 1e3 for x in XC], ux, lw=1.0,
                         label=f"t = {t*1e6:.2f} $\\mu$s")
            ax2.set_xlabel("x [mm]")
            ax2.set_ylabel("$u_x$ [m/s]")
            ax2.set_title(f"{mode}: $u_p$ = {r['up_nominal']:.0f} m/s "
                          "wave profiles")
            ax2.legend(fontsize=7)
            ax2.grid(alpha=0.3)
            fig2.tight_layout()
            fig2.savefig(f"profiles_{mode}.png", dpi=180)
            print(f"wrote profiles_{mode}.png")
    except ImportError:
        print("\n(matplotlib not available: skipping figures)")

    if passed is False:
        sys.exit(1)


if __name__ == "__main__":
    main()
