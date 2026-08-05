#!/usr/bin/env python3
# =========================================================================
#  validate_riemann.py
#
#  Solid Riemann (Sod-type) problem for crystaldynamicsFoam's Mie-Gruneisen EOS mode.
#
#  Initial data (both halves AT REST, diaphragm at x0 = 5 mm):
#    left  : the principal-Hugoniot state of up_L = 800 m/s, at rest --
#            rho_L = rho0(1+mu_L), e_L = up_L^2/2, p_L = 36.66 GPa, HOT
#    right : ambient (rho0, e = 0, p = 0)
#
#  Exact solution: a right-running SHOCK into the ambient material, a
#  left-running RELEASE FAN into the hot compressed material, and between
#  them a CONTACT where p and u are continuous but rho and T jump --
#  because the two sides reach the same (p*, u*) along different
#  thermodynamic paths (single shock vs shock-then-release).  That
#  entropy contact is the observable that uniquely certifies the energy
#  bookkeeping: neither the Hugoniot suite nor the release benchmark
#  tests T directly.
#
#  Star state from impedance matching:
#      right shock locus : p_s(u) = rho0 (c0 + s u) u
#      left release locus: p_r(u) from the isentrope through the left
#          state, with du/dmu = -c/(1+mu) (left-facing fan: u rises as
#          the material expands toward the low-pressure side)
#      p_r(u*) = p_s(u*)  ->  (u*, p*)
#
#  The legacy closure (p = K mu, c = const) is solved with the same
#  structure for comparison; note it predicts NO density jump at the
#  contact (same p -> same mu -> same rho on both sides).
#
#  Temperature convention: the exact T uses the SAME leading-order
#  recovery as the solver, T = Tref + (e - e_H(mu))/Cv, so this is
#  verification against the model's own exact solution.
#
#  Pure standard library; matplotlib used for figures if importable.
# =========================================================================

import os
import re
import math
import sys

# ----------------------- case constants (must match caseTemplate) --------
L_BAR   = 1.0e-2
N_CELLS = 1000
DX      = L_BAR / N_CELLS
XC      = [(i + 0.5) * DX for i in range(N_CELLS)]
X0      = 5.0e-3          # diaphragm position


def _set_mesh(n_cells, l_bar):
    """Rebind the mesh globals (cell count / bar length / centres)."""
    global N_CELLS, L_BAR, DX, XC
    if n_cells == N_CELLS and abs(l_bar - L_BAR) < 1e-12:
        return
    N_CELLS, L_BAR = n_cells, l_bar
    DX = L_BAR / N_CELLS
    XC = [(i + 0.5) * DX for i in range(N_CELLS)]
    print(f"  (mesh configured: {N_CELLS} cells, L = {L_BAR*1e3:.2f} mm, "
          f"dx = {DX*1e6:.2f} um)")


def configure_mesh(case):
    """Read L and the x cell count from the run's own blockMeshDict, so
    the validator follows any resolution/length change automatically.
    Falls back to the current values (plus the self-healing length check
    in read_internal_field) if parsing fails."""
    bmd = os.path.join(case, "system", "blockMeshDict")
    if not os.path.isfile(bmd):
        return
    with open(bmd) as fh:
        txt = fh.read()
    txt = re.sub(r"/\*.*?\*/", " ", txt, flags=re.S)
    txt = re.sub(r"//[^\n]*", " ", txt)
    mL = re.search(r"\bL\s+([0-9.eE+\-]+)\s*;", txt)
    mN = re.search(r"hex\s*\([^)]*\)\s*\(\s*(\d+)\s+\d+\s+\d+\s*\)", txt)
    if mN:
        _set_mesh(int(mN.group(1)),
                  float(mL.group(1)) if mL else L_BAR)

RHO0    = 8930.0
C0      = 3940.0
S       = 1.489
G0      = 1.99
CV      = 385.0
TREF    = 300.0
C11, C12 = 168.4e9, 121.4e9
K_BULK  = (C11 + 2.0 * C12) / 3.0
CK      = math.sqrt(K_BULK / RHO0)

UP_L    = 800.0           # Hugoniot parameterisation of the LEFT state

TOL_U   = 0.03            # PASS tolerances (MieGruneisen)
TOL_P   = 0.04
TOL_US  = 0.03
TOL_DT  = 0.20            # contact temperature jump
TOL_CONT_U = 0.015        # contact continuity, fraction of u*
TOL_CONT_P = 0.02         # contact continuity, fraction of p*

RUNS_DIR = "runs"

# ----------------------- Mie-Gruneisen algebra ---------------------------

A0 = RHO0 * C0 * C0

def pH(mu):
    if mu >= 0.0:
        D = 1.0 - (S - 1.0) * mu
        return A0 * mu * (1.0 + mu) / (D * D)
    return A0 * mu

def eH(mu):
    if mu >= 0.0:
        return 0.5 * pH(mu) * mu / (RHO0 * (1.0 + mu))
    return 0.5 * (A0 / RHO0) * mu * mu

def dpH(mu):
    if mu >= 0.0:
        D = 1.0 - (S - 1.0) * mu
        return A0 * ((1.0 + 2.0 * mu) * D
                     + 2.0 * (S - 1.0) * mu * (1.0 + mu)) / (D ** 3)
    return A0

def deH(mu, h=1e-7):
    return (eH(mu + h) - eH(mu - h)) / (2.0 * h)

def p_of(mu, e):
    return pH(mu) + G0 * RHO0 * (e - eH(mu))

def c_isen(mu, e):
    dpdmu_s = dpH(mu) - G0 * RHO0 * deH(mu) \
              + G0 * p_of(mu, e) / (1.0 + mu) ** 2
    return math.sqrt(max(dpdmu_s, 1.0) / RHO0)

def T_of(mu, e):
    """Leading-order recovery -- SAME convention as the solver."""
    return max(TREF + (e - eH(max(mu, 0.0))) / CV, 1.0) if mu >= 0.0 \
        else max(TREF + (e - 0.5 * (A0 / RHO0) * mu * mu) / CV, 1.0)

def hugoniot_state(up):
    Us = C0 + S * up
    mu = up / (Us - up)
    return Us, mu, RHO0 * Us * up, 0.5 * up * up

def left_release_path(mu_L, e_L, nstep=6000):
    """Isentrope through the left state, LEFT-facing fan: u rises as mu
    falls.  Returns samples [(u, p, c, mu, e, T)], u ascending from 0."""
    dmu = mu_L / nstep

    def rhs(mu, e):
        p = p_of(mu, e)
        return p / (RHO0 * (1.0 + mu) ** 2), c_isen(mu, e) / (1.0 + mu)

    mu, e, u = mu_L, e_L, 0.0
    path = [(u, p_of(mu, e), c_isen(mu, e), mu, e, T_of(mu, e))]
    while mu > dmu:
        k1e, k1u = rhs(mu, e)
        k2e, k2u = rhs(mu - 0.5 * dmu, e - 0.5 * dmu * k1e)
        k3e, k3u = rhs(mu - 0.5 * dmu, e - 0.5 * dmu * k2e)
        k4e, k4u = rhs(mu - dmu,       e - dmu * k3e)
        e  -= dmu * (k1e + 2 * k2e + 2 * k3e + k4e) / 6.0
        u  += dmu * (k1u + 2 * k2u + 2 * k3u + k4u) / 6.0   # u RISES
        mu -= dmu
        path.append((u, p_of(mu, e), c_isen(mu, e), mu, e, T_of(mu, e)))
    return path

def path_at_u(path, u_t):
    for a, b in zip(path, path[1:]):
        if (a[0] - u_t) * (b[0] - u_t) <= 0.0 and a[0] != b[0]:
            f = (u_t - a[0]) / (b[0] - a[0])
            return tuple(ai + f * (bi - ai) for ai, bi in zip(a, b))
    return path[-1]

def solve_star_mg():
    UsL_shk, mu_L, p_L, e_L = hugoniot_state(UP_L)
    path = left_release_path(mu_L, e_L)

    def mismatch(u):
        p_rel = path_at_u(path, u)[1]
        p_shk = RHO0 * (C0 + S * u) * u
        return p_rel - p_shk

    lo, hi = 1.0, UP_L - 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mismatch(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    u_star = 0.5 * (lo + hi)
    st = path_at_u(path, u_star)          # left-of-contact star state
    p_star = st[1]
    Us = C0 + S * u_star
    muR = u_star / (Us - u_star)
    return dict(
        mu_L=mu_L, p_L=p_L, e_L=e_L, c_L=path[0][2], T_L=path[0][5],
        path=path, u=u_star, p=p_star, Us=Us,
        rhoR=RHO0 * (1.0 + muR), TR=T_of(muR, eH(muR)),
        rhoL=RHO0 * (1.0 + st[3]), TL=st[5], cLs=st[2],
    )

# ----------------------- legacy closure exact -----------------------------

def solve_star_legacy():
    """p = K mu, c = CK constant.  Left state: p_L = K mu_L at rest."""
    _, mu_L, _, _ = hugoniot_state(UP_L)   # same rho_L initialisation
    p_L = K_BULK * mu_L

    def u_rel(mu):
        return CK * math.log((1.0 + mu_L) / (1.0 + mu))

    def mismatch(mu):
        return u_rel(mu) - CK * mu / math.sqrt(1.0 + mu)

    lo, hi = 1e-6, mu_L
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mismatch(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    mu_star = 0.5 * (lo + hi)
    u_star = CK * mu_star / math.sqrt(1.0 + mu_star)
    return dict(
        mu_L=mu_L, p_L=p_L, c_L=CK, u=u_star, p=K_BULK * mu_star,
        Us=CK * math.sqrt(1.0 + mu_star),
        rhoR=RHO0 * (1.0 + mu_star), rhoL=RHO0 * (1.0 + mu_star),  # NO jump
    )

# ----------------------- exact profiles -----------------------------------

def exact_profile_mg(sol, t):
    """(u, p, rho, T)(x) at time t for the MG solution."""
    path = sol["path"]
    head = -sol["c_L"]                       # lab speed of fan head
    tail = sol["u"] - sol["cLs"]             # lab speed of fan tail
    xc_contact = X0 + sol["u"] * t
    xc_shock = X0 + sol["Us"] * t

    def fan_state(xi):
        lo, hi = 0.0, sol["u"]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            st = path_at_u(path, mid)
            if mid - st[2] < xi:
                lo = mid
            else:
                hi = mid
        return path_at_u(path, 0.5 * (lo + hi))

    u, p, rho, T = [], [], [], []
    for x in XC:
        xi = (x - X0) / t
        if xi <= head:
            st = (0.0, sol["p_L"], 0.0, sol["mu_L"], 0.0, sol["T_L"])
        elif xi <= tail:
            st = fan_state(xi)
        elif x <= xc_contact:
            st = path_at_u(path, sol["u"])
        elif x <= xc_shock:
            st = (sol["u"], sol["p"], 0.0,
                  sol["rhoR"] / RHO0 - 1.0, 0.0, sol["TR"])
        else:
            st = (0.0, 0.0, 0.0, 0.0, 0.0, TREF)
        u.append(st[0])
        p.append(st[1])
        rho.append(RHO0 * (1.0 + st[3]))
        T.append(st[5])
    return u, p, rho, T

def exact_profile_legacy(sol, t):
    mu_L = sol["mu_L"]
    head, tail = -CK, sol["u"] - CK
    xc_contact, xc_shock = X0 + sol["u"] * t, X0 + sol["Us"] * t
    u, p, rho, T = [], [], [], []
    for x in XC:
        xi = (x - X0) / t
        if xi <= head:
            uu, mm = 0.0, mu_L
        elif xi <= tail:
            uu = xi + CK
            mm = (1.0 + mu_L) / math.exp(uu / CK) - 1.0
        elif x <= xc_shock:
            uu, mm = sol["u"], sol["rhoR"] / RHO0 - 1.0
        else:
            uu, mm = 0.0, 0.0
        u.append(uu)
        p.append(K_BULK * mm if x <= xc_shock else 0.0)
        rho.append(RHO0 * (1.0 + mm))
        # exact legacy T: total-energy balance with the athermal closure.
        # Release extracts the expansion work from Cv*T (work integral has
        # the closed form below -> NEGATIVE T in the released region);
        # the shock books ALL Hugoniot energy as heat.
        def _W(m):
            return (K_BULK / RHO0) * (math.log(1.0 + m) + 1.0 / (1.0 + m))
        if xi <= head:
            T.append(TREF)
        elif x <= xc_contact:
            T.append(TREF - (_W(mu_L) - _W(mm)) / CV)
        elif x <= xc_shock:
            T.append(TREF + 0.5 * sol["u"] ** 2 / CV)
        else:
            T.append(TREF)
    return u, p, rho, T

# ----------------------- OpenFOAM ASCII field parsing ---------------------

def _matched_block(text, start):
    depth = 0
    for i in range(start, len(text)):
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
    raise ValueError("unbalanced parentheses")

def read_internal_field(fname, ncomp):
    with open(fname) as fh:
        txt = fh.read()
    sub = txt[txt.find("internalField"):]
    m = re.match(r"internalField\s+(uniform|nonuniform)", sub)
    if m.group(1) == "uniform":
        if ncomp == 1:
            v = re.match(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", sub)
            return [float(v.group(1))] * N_CELLS
        v = re.match(r"internalField\s+uniform\s+\(([^)]*)\)\s*;", sub)
        comps = [float(t) for t in v.group(1).split()]
        return [comps[:] for _ in range(N_CELLS)]
    mn = re.match(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*", sub)
    n_file = int(mn.group(1))
    if n_file != N_CELLS:
        # resolution changed and blockMeshDict parsing did not catch it:
        # trust the field file (bar length assumed unchanged)
        _set_mesh(n_file, L_BAR)
    body = _matched_block(sub, sub.find("(", mn.end()))
    if ncomp == 1:
        return [float(t) for t in body.split()]
    return [[float(t) for t in tup.split()]
            for tup in re.findall(r"\(([^)]+)\)", body)]

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

def window_mean(arr, xlo, xhi):
    vals = [arr[i] for i, x in enumerate(XC) if xlo <= x <= xhi]
    return sum(vals) / len(vals) if vals else float("nan")

def front_fit(case, threshold, t_min_frac=0.3):
    pts = []
    for t, d in time_dirs(case):
        uf = os.path.join(d, "U")
        if not os.path.isfile(uf):
            continue
        ux = [v[0] for v in read_internal_field(uf, 3)]
        idx = [i for i, v in enumerate(ux) if v >= threshold and XC[i] > X0]
        if idx and max(idx) < N_CELLS - 1:
            i0 = max(idx)
            f = (ux[i0] - threshold) / (ux[i0] - ux[i0 + 1]) \
                if ux[i0] != ux[i0 + 1] else 0.0
            pts.append((t, XC[i0] + f * DX))
    if len(pts) < 4:
        return None
    t_last = pts[-1][0]
    sel = [(t, x) for t, x in pts if t >= t_min_frac * t_last]
    n = len(sel)
    sx = sum(t for t, x in sel)
    sy = sum(x for t, x in sel)
    sxx = sum(t * t for t, x in sel)
    sxy = sum(t * x for t, x in sel)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)

# ----------------------- main ---------------------------------------------

def main():
    mg = solve_star_mg()
    lg = solve_star_legacy()

    print("=" * 72)
    print("Exact Riemann solutions (left = Hugoniot(800 m/s) at rest | ambient)")
    print("=" * 72)
    print(f"  left state    : rho_L = {RHO0*(1+mg['mu_L']):.1f} kg/m3, "
          f"e_L = {mg['e_L']:.0f} J/kg")
    print(f"                  p_L(MG) = {mg['p_L']/1e9:.2f} GPa | "
          f"p_L(legacy) = {lg['p_L']/1e9:.2f} GPa")
    print(f"  MieGruneisen  : u* = {mg['u']:6.1f} m/s   p* = "
          f"{mg['p']/1e9:6.2f} GPa   Us = {mg['Us']:6.1f} m/s")
    print(f"                  contact: rho {mg['rhoL']:.0f} | "
          f"{mg['rhoR']:.0f} kg/m3 (jump {mg['rhoL']-mg['rhoR']:+.0f}),  "
          f"T {mg['TL']:.0f} | {mg['TR']:.0f} K "
          f"(jump {mg['TL']-mg['TR']:+.0f} K)")
    print(f"  linearElastic : u* = {lg['u']:6.1f} m/s   p* = "
          f"{lg['p']/1e9:6.2f} GPa   Us = {lg['Us']:6.1f} m/s   "
          "(NO contact density jump)")
    print(f"  discrimination: p* {100*(mg['p']/lg['p']-1):+.0f}%,  "
          f"u* {100*(mg['u']/lg['u']-1):+.0f}%")

    if not os.path.isdir(RUNS_DIR):
        sys.exit("\nNo runs/ directory found -- run ./Allrun first.")

    results = {}
    for mode in sorted(os.listdir(RUNS_DIR)):
        case = os.path.join(RUNS_DIR, mode, "riemann")
        if not os.path.isdir(case):
            continue
        configure_mesh(case)
        tds = time_dirs(case)
        if not tds:
            print(f"\n{mode}: no output; skipped")
            continue
        t_fin, d_fin = tds[-1]

        ux = [v[0] for v in read_internal_field(os.path.join(d_fin, "U"), 3)]
        pp = read_internal_field(os.path.join(d_fin, "p"), 1)
        rr = read_internal_field(os.path.join(d_fin, "rho"), 1)
        TT = read_internal_field(os.path.join(d_fin, "T"), 1) \
            if os.path.isfile(os.path.join(d_fin, "T")) else [float("nan")] * N_CELLS

        sol = mg if mode == "MieGruneisen" else lg
        x_tail = X0 + (sol["u"] - (sol["cLs"] if mode == "MieGruneisen"
                                   else CK)) * t_fin
        x_cont = X0 + sol["u"] * t_fin
        x_shck = X0 + sol["Us"] * t_fin
        lw = (x_tail + 0.8e-3, x_cont - 0.4e-3)
        rw = (x_cont + 0.4e-3, x_shck - 0.4e-3)

        meas = dict(
            t=t_fin, ux=ux, p=pp, rho=rr, T=TT,
            uL=window_mean(ux, *lw), uR=window_mean(ux, *rw),
            pL=window_mean(pp, *lw), pR=window_mean(pp, *rw),
            rhoL=window_mean(rr, *lw), rhoR=window_mean(rr, *rw),
            TL=window_mean(TT, *lw), TR=window_mean(TT, *rw),
            Us=front_fit(case, 0.5 * sol["u"]),
        )
        results[mode] = meas

        print(f"\n=== mode: {mode}  (final time {t_fin*1e6:.2f} us) ===")
        print(f"  star plateaus : u = {meas['uL']:6.1f} | {meas['uR']:6.1f} "
              f"m/s (exact {sol['u']:.1f}),  p = {meas['pL']/1e9:5.2f} | "
              f"{meas['pR']/1e9:5.2f} GPa (exact {sol['p']/1e9:.2f})")
        print(f"  shock speed   : {meas['Us']:6.1f} m/s "
              f"(exact {sol['Us']:.1f})" if meas['Us'] else
              "  shock speed   : NOT MEASURED")
        print(f"  contact       : drho = {meas['rhoL']-meas['rhoR']:+6.1f} "
              f"kg/m3, dT = {meas['TL']-meas['TR']:+6.1f} K")

        if mode == "MieGruneisen":
            u_m = 0.5 * (meas["uL"] + meas["uR"])
            p_m = 0.5 * (meas["pL"] + meas["pR"])
            e_u = abs(u_m - sol["u"]) / sol["u"]
            e_p = abs(p_m - sol["p"]) / sol["p"]
            e_us = abs(meas["Us"] - sol["Us"]) / sol["Us"] \
                if meas["Us"] else 1.0
            dT_ex = sol["TL"] - sol["TR"]
            e_dt = abs((meas["TL"] - meas["TR"]) - dT_ex) / dT_ex
            c_u = abs(meas["uL"] - meas["uR"]) / sol["u"]
            c_p = abs(meas["pL"] - meas["pR"]) / sol["p"]
            ok = (e_u < TOL_U and e_p < TOL_P and e_us < TOL_US
                  and e_dt < TOL_DT and c_u < TOL_CONT_U
                  and c_p < TOL_CONT_P)
            print(f"  -> errors: u* {100*e_u:.1f}% (tol {100*TOL_U:.0f}), "
                  f"p* {100*e_p:.1f}% (tol {100*TOL_P:.0f}), "
                  f"Us {100*e_us:.1f}% (tol {100*TOL_US:.0f}), "
                  f"contact dT {100*e_dt:.0f}% (tol {100*TOL_DT:.0f})")
            print(f"  -> contact continuity: |du| {100*c_u:.2f}% of u*, "
                  f"|dp| {100*c_p:.2f}% of p*")
            print(f"  -> {'PASS' if ok else 'FAIL'}: EOS mode "
                  f"{'matches' if ok else 'does NOT match'} the exact "
                  "Riemann solution (incl. the entropy contact)")
            results[mode]["pass"] = ok
        else:
            print("  -> vs the REAL (MG) solution: star pressure "
                  f"{100*(meas['pR']/mg['p']-1):+.0f}%, star velocity "
                  f"{100*(meas['uR']/mg['u']-1):+.0f}% -- and no physical "
                  "contact density/temperature jump: the missing physics")

    # ------------------------- figures -----------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for mode, meas in results.items():
            sol = mg if mode == "MieGruneisen" else lg
            ex = exact_profile_mg(sol, meas["t"]) if mode == "MieGruneisen" \
                else exact_profile_legacy(sol, meas["t"])
            fig, axs = plt.subplots(2, 2, figsize=(9.6, 6.6), sharex=True)
            xs = [x * 1e3 for x in XC]
            panels = [("$u_x$ [m/s]", meas["ux"], ex[0], 1.0),
                      ("p [GPa]", meas["p"], ex[1], 1e-9),
                      ("$\\rho$ [kg/m$^3$]", meas["rho"], ex[2], 1.0),
                      ("T [K]", meas["T"], ex[3], 1.0)]
            for ax, (lab, sim, exact, sc) in zip(axs.flat, panels):
                ax.plot(xs, [v * sc for v in exact], "k-", lw=1.2,
                        label="exact")
                ax.plot(xs, [v * sc for v in sim], ".", ms=2.5,
                        color="tab:red" if mode == "MieGruneisen"
                        else "tab:blue", label="crystaldynamicsFoam")
                ax.set_ylabel(lab)
                ax.grid(alpha=0.3)
            axs[0, 0].legend(fontsize=8)
            for ax in axs[1]:
                ax.set_xlabel("x [mm]")
            fig.suptitle(f"Solid Riemann problem, {mode}, "
                         f"t = {meas['t']*1e6:.2f} $\\mu$s")
            fig.tight_layout()
            fig.savefig(f"riemann_profiles_{mode}.png", dpi=180)
            print(f"\nwrote riemann_profiles_{mode}.png")

        # p-u diagram
        fig, ax = plt.subplots(figsize=(6.6, 4.8))
        uu = [i * mg["u"] * 1.6 / 200 for i in range(201)]
        ax.plot(uu, [RHO0 * (C0 + S * u) * u / 1e9 for u in uu], "k-",
                label="right-shock Hugoniot locus (MG)")
        ax.plot([st[0] for st in mg["path"]],
                [st[1] / 1e9 for st in mg["path"]], "k-.",
                label="left release isentrope (MG)")
        ax.plot(mg["u"], mg["p"] / 1e9, "r*", ms=14,
                label=f"MG star ({mg['u']:.0f} m/s, {mg['p']/1e9:.1f} GPa)")
        ax.plot(lg["u"], lg["p"] / 1e9, "b*", ms=14,
                label=f"legacy star ({lg['u']:.0f} m/s, {lg['p']/1e9:.1f} GPa)")
        for mode, meas in results.items():
            col = "tab:red" if mode == "MieGruneisen" else "tab:blue"
            ax.plot(0.5 * (meas["uL"] + meas["uR"]),
                    0.5 * (meas["pL"] + meas["pR"]) / 1e9, "o", ms=9,
                    mfc="none", mew=1.8, color=col,
                    label=f"{mode} measured")
        ax.set_xlabel("particle velocity u [m/s]")
        ax.set_ylabel("pressure p [GPa]")
        ax.set_title("Riemann p-u diagram: loci intersection = star state")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig("riemann_pu_diagram.png", dpi=180)
        print("wrote riemann_pu_diagram.png")
    except ImportError:
        print("\n(matplotlib not available: skipping figures)")

    if results.get("MieGruneisen", {}).get("pass") is False:
        sys.exit(1)


if __name__ == "__main__":
    main()
