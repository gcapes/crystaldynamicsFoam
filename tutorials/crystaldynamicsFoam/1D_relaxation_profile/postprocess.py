#!/usr/bin/env python3
"""
Post-processing for equilibrium double-obstacle interface profile case.

Produces:
  profile_evolution.png      -- one panel per mesh resolution (2/4/6/8/10/20
                                 cells across eta), phi_1(x) at every
                                 timestep vs analytical sinusoid (opacity =
                                 how close to converged)
  profile_and_error_vs_h.png -- (a) converged profile, coarsest vs finest
                                 mesh, vs analytical; (b) L2/Linf/energy
                                 error vs mesh resolution h = eta/dx
  convergence_table.txt      -- L2/Linf profile error, energy-integral
                                 relative error, and the dt-halving /
                                 ramp-IC cross-checks

Run after Allrun (which calls this automatically), or standalone from this
directory: `python3 postprocess.py`.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.collections import LineCollection

import foam_io

HERE = os.path.dirname(os.path.abspath(__file__))

MESH_CASES = [("mesh2", 2), ("mesh4", 4), ("mesh6", 6),
              ("mesh8", 8), ("mesh10", 10), ("mesh20", 20)]

PANEL_LABELS = "abcdefghij"


def analytic_profile(x, eta):
    x = np.asarray(x)
    return np.where(x <= -eta / 2, 1.0, np.where(x >= eta / 2, 0.0,
                    0.5 * (1.0 - np.sin(np.pi * x / eta))))


def _antiderivative(u, eta):
    """Antiderivative of analytic_profile, continuous across both kinks."""
    u = np.asarray(u)
    left = u
    mid = 0.5 * u + (eta / (2 * np.pi)) * np.cos(np.pi * u / eta) - eta / 4.0
    right = np.zeros_like(u)
    return np.where(u <= -eta / 2, left, np.where(u >= eta / 2, right, mid))


def analytic_cell_average(xc, dx, eta):
    """Exact cell-average of the analytical profile over [x-dx/2, x+dx/2] --
    the correct comparison for a finite-volume field: a cell-centre point
    sample is sensitive to whether the compact-support kink at x=+-eta/2
    lands on a cell centre or a face, which can make the pointwise error
    non-monotonic under refinement for reasons unrelated to accuracy."""
    xc = np.asarray(xc)
    return (_antiderivative(xc + dx / 2, eta) - _antiderivative(xc - dx / 2, eta)) / dx


def interfacial_energy_density(phi, x, gamma, eta):
    """f(x) = (4*gamma/eta) * [ (eta^2/pi^2)*(dphi/dx)^2 + phi*(1-phi) ]"""
    dphidx = np.gradient(phi, x)
    return (4.0 * gamma / eta) * ((eta ** 2 / np.pi ** 2) * dphidx ** 2 + phi * (1 - phi))


def total_interfacial_energy(phi, x, gamma, eta):
    f = interfacial_energy_density(phi, x, gamma, eta)
    return float(np.trapz(f, x))


def converged_profile(case_dir):
    params = foam_io.load_params(case_dir)
    times = foam_io.list_times(case_dir)
    assert times, f"no time directories found in {case_dir}"
    t_final, name_final = times[-1]
    xc = np.array(foam_io.cell_centres(params))
    phi = np.array(foam_io.expand(
        foam_io.read_scalar_field(os.path.join(case_dir, name_final, "n.0")),
        params["Nx"]))
    return params, xc, phi, t_final


def l2_linf_error(phi, phi_ana):
    """L2 = RMS(phi_numeric - phi_analytical) over all cells; Linf = largest
    single-cell difference."""
    diff = phi - phi_ana
    l2 = float(np.sqrt(np.mean(diff ** 2)))
    linf = float(np.max(np.abs(diff)))
    return l2, linf


def plot_profile_evolution(mesh_cases, colors):
    ncols = 3
    nrows = math.ceil(len(mesh_cases) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.6 * nrows),
                              sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    alpha_min = 0.06

    tconv_over_tau = {}
    for panel_idx, ((case_dir, n), color, ax) in enumerate(zip(mesh_cases, colors, axes)):
        params = foam_io.load_params(case_dir)
        eta, dx = params["eta"], params["dx"]
        end_time, tau = params["end_time"], params["tau"]
        xc = np.array(foam_io.cell_centres(params)) / eta
        times = foam_io.list_times(case_dir)
        n_t = len(times)

        phis = [np.array(foam_io.expand(
            foam_io.read_scalar_field(os.path.join(case_dir, name, "n.0")),
            params["Nx"])) for _, name in times]

        segments = []
        rgba = np.tile(np.array(color), (n_t, 1))
        for i, ((t, name), phi) in enumerate(zip(times, phis)):
            segments.append(np.column_stack([xc, phi]))
            # opacity tracks simulated time (not snapshot index), so panels
            # -- all sharing the same end_time = 3*tau -- are comparable
            frac = np.clip(t / end_time, 0.0, 1.0)
            rgba[i, 3] = alpha_min + (1.0 - alpha_min) * frac

        # convergence time: first snapshot after which ||phi(t)-phi(final)||_2
        # stays below tolerance for the rest of the run
        phis_arr = np.array(phis)
        l2_to_final = np.sqrt(np.mean((phis_arr - phis_arr[-1]) ** 2, axis=1))
        conv_tol = 1e-3
        idx_conv = next((i for i in range(n_t) if np.all(l2_to_final[i:] < conv_tol)), n_t - 1)
        tconv_over_tau[n] = times[idx_conv][0] / tau

        # kept fully vector (not rasterized): matplotlib 3.6.3 partially drops
        # per-segment-alpha rasterized LineCollections in PDF output
        lc = LineCollection(segments, colors=rgba, linewidths=1.0, zorder=1)
        ax.add_collection(lc)

        # converged profile as hollow dots so the analytical line (on top,
        # drawn next) stays visible through their centres
        ax.plot(xc, np.array(foam_io.expand(
            foam_io.read_scalar_field(os.path.join(case_dir, times[-1][1], "n.0")),
            params["Nx"])), "o", mfc="none", mec=color, mew=1.2, ms=4.5, zorder=2)

        xs_fine = np.linspace(-2.0, 2.0, 2000)
        ax.plot(xs_fine, analytic_profile(xs_fine, 1.0),
                 "k-", lw=1.1, alpha=0.85, zorder=3, label="analytical")

        ax.set_xlim(-2, 2)
        ax.set_ylim(-0.05, 1.05)
        ax.text(-0.08, 1.12, f"{PANEL_LABELS[panel_idx]})", transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="left")
        ax.grid(alpha=0.3)
        ax.set_title(f"{n} cells/$\\eta$ (dx = {dx:.2g} m)", fontstyle="italic", fontsize=10)

    # report convergence time per panel if it varies noticeably mesh-to-mesh
    tconv_vals = list(tconv_over_tau.values())
    if (max(tconv_vals) - min(tconv_vals)) > 0.05 * np.mean(tconv_vals):
        for (case_dir, n), ax in zip(mesh_cases, axes):
            base = ax.get_title()
            ax.set_title(base[:-1] + f", $t_{{conv}}\\approx${tconv_over_tau[n]:.2f}$\\tau$)",
                         fontstyle="italic", fontsize=10)

    for ax in axes[len(mesh_cases):]:
        ax.axis("off")
    for ax in axes[:len(mesh_cases)]:
        if ax.get_subplotspec().is_last_row():
            ax.set_xlabel(r"$x/\eta$")
        if ax.get_subplotspec().is_first_col():
            ax.set_ylabel(r"$\phi_1$")

    handles = [plt.Line2D([], [], color="k", ls="-", lw=1.1, alpha=0.85, label="analytical")]
    fig.legend(handles=handles, loc="lower center", ncol=1, fontsize=9,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(os.path.join(HERE, "profile_evolution.png"), dpi=150,
                bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    return tconv_over_tau


def plot_profile_and_error_vs_h(mesh_cases, colors, case_data):
    """Two-panel summary: (a) converged profile, coarsest vs finest mesh,
    against the analytical sinusoid; (b) L2/Linf/energy error vs h."""
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.0, 4.0))

    ns = [n for _, n in mesh_cases]
    h_lo, h_hi = min(ns), max(ns)
    markers = {h_lo: "o", h_hi: "x"}
    zorder_by_h = {h_hi: 2, h_lo: 3}

    xs_fine = np.linspace(-1.5, 1.5, 2000)
    axa.plot(xs_fine, analytic_profile(xs_fine, 1.0), "k-", lw=1.3, zorder=1,
              label="analytical")
    for h, color in zip(ns, colors):
        if h not in (h_lo, h_hi):
            continue
        xc, phi, eta = case_data[h]["xc"], case_data[h]["phi"], case_data[h]["eta"]
        axa.plot(xc / eta, phi, markers[h], ls="none", mfc="none", mec=color,
                  mew=1.1, ms=3.5, label=f"$h = {h}$", zorder=zorder_by_h[h])

    axa.set_xlim(-1.5, 1.5)
    axa.set_ylim(-0.05, 1.05)
    axa.set_xlabel(r"$x/\eta$")
    axa.set_ylabel(r"$\varphi_1$")
    axa.grid(alpha=0.3)
    axa.legend(fontsize=9, loc="upper right", framealpha=0.9)
    axa.text(-0.12, 1.08, "(a)", transform=axa.transAxes, fontsize=12,
              fontweight="bold", va="top", ha="left")

    ns_arr = np.array(ns)
    l2s_pct = np.array([case_data[n]["l2"] for n in ns]) * 100
    linfs_pct = np.array([case_data[n]["linf"] for n in ns]) * 100
    energy_errs_pct = np.array([case_data[n]["energy_err"] for n in ns]) * 100

    axb.semilogy(ns_arr, l2s_pct, "o-",
                 label=r"$\sqrt{\langle(\varphi_1-\varphi_1^{\mathrm{ana}})^2\rangle}$")
    axb.semilogy(ns_arr, linfs_pct, "s-",
                 label=r"$\max|\varphi_1-\varphi_1^{\mathrm{ana}}|$")
    axb.semilogy(ns_arr, energy_errs_pct, "^-",
                 label=r"$\left|\dfrac{\int f\,dx-\gamma}{\gamma}\right|$")

    for tol in (5.0, 1.0):
        axb.axhline(tol, color="0.7", lw=1.0, ls="--", zorder=0)
        axb.text(1.02, tol, f"{tol:g}%", color="0.5", fontsize=8, va="center",
                  transform=axb.get_yaxis_transform())

    axb.set_xlabel(r"$h = \eta/\Delta x$")
    axb.set_ylabel("error (%)")
    axb.set_xticks(ns_arr)
    axb.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f"{y:g}%"))
    axb.legend(fontsize=8.5, loc="upper right", framealpha=0.9)
    axb.grid(alpha=0.3, which="both")
    axb.text(-0.12, 1.08, "(b)", transform=axb.transAxes, fontsize=12,
              fontweight="bold", va="top", ha="left")

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "profile_and_error_vs_h.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)


def main():
    mesh_cases = [(os.path.join(HERE, name), n) for name, n in MESH_CASES
                  if os.path.isdir(os.path.join(HERE, name))]
    if not mesh_cases:
        print("No mesh2/4/6/8/10/20 case directories found -- run Allrun first.")
        return

    report_lines = []

    def log(line=""):
        print(line)
        report_lines.append(line)

    log("Case 1 -- equilibrium double-obstacle interface profile")
    log("=========================================================")

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(mesh_cases)))

    tconv_over_tau = plot_profile_evolution(mesh_cases, colors)
    log("Saved profile_evolution.png")

    tconv_vals = list(tconv_over_tau.values())
    if (max(tconv_vals) - min(tconv_vals)) > 0.05 * np.mean(tconv_vals):
        log(f"\nConvergence time (||phi(t)-phi(final)||_2 < {1e-3:.0e}, held to the end) "
            f"varies by mesh -- shown per panel:")
        for n, tc in tconv_over_tau.items():
            log(f"  n={n:3d}: t_conv/tau = {tc:.4f}")
    else:
        log(f"\nConvergence time is essentially the same across meshes: "
            f"t_conv/tau ~ {np.mean(tconv_vals):.2f}")

    # Converged-profile comparison against the cell-averaged analytical
    # profile (see analytic_cell_average) -- also builds case_data, reused
    # below for the profile_and_error_vs_h.png panels.
    log("\nConverged-profile comparison against analytical sinusoid")
    log("  L2   = root-mean-square(phi_numeric - phi_analytical) over all cells")
    log("  Linf = largest single-cell |phi_numeric - phi_analytical|")
    log("  energy rel. err = (integrated interfacial energy - gamma) / gamma")
    log(f"{'mesh (cells/eta)':>18} {'dx/eta':>10} {'t/tau':>8} "
        f"{'L2 error':>12} {'Linf error':>12} {'energy rel. err':>16}")

    case_data = {}
    for case_dir, n in mesh_cases:
        params, xc, phi, t_final = converged_profile(case_dir)
        eta, gamma, tau, dx = params["eta"], params["gamma"], params["tau"], params["dx"]

        phi_ana_avg = analytic_cell_average(xc, dx, eta)
        l2, linf = l2_linf_error(phi, phi_ana_avg)
        F = total_interfacial_energy(phi, xc, gamma, eta)
        energy_rel_err = (F - gamma) / gamma

        log(f"{n:>18d} {1.0/n:>10.4f} {t_final/tau:>8.2f} "
            f"{l2:>12.4e} {linf:>12.4e} {energy_rel_err:>16.4e}")

        case_data[n] = dict(xc=xc, phi=phi, eta=eta, l2=l2, linf=linf,
                             energy_err=abs(energy_rel_err))

    log("\nCoarsest resolution meeting a given energy-error tolerance:")
    for tol in (0.05, 0.01):
        ok = [n for n, d in case_data.items() if d["energy_err"] < tol]
        verdict = f"n={min(ok)} cells/eta" if ok else "none of the tested resolutions"
        log(f"  <{tol:.0%}: {verdict}")

    plot_profile_and_error_vs_h(mesh_cases, colors, case_data)
    log("Saved profile_and_error_vs_h.png")

    # Monotonic-decrease check on total interfacial free energy F(t)
    log("\nMonotonic-decrease check on total interfacial free energy F(t)")
    for case_dir, n in mesh_cases:
        params = foam_io.load_params(case_dir)
        eta, gamma, tau = params["eta"], params["gamma"], params["tau"]
        xc = np.array(foam_io.cell_centres(params))
        times = foam_io.list_times(case_dir)

        Fs = []
        for t, name in times:
            phi = np.array(foam_io.expand(
                foam_io.read_scalar_field(os.path.join(case_dir, name, "n.0")),
                params["Nx"]))
            Fs.append(total_interfacial_energy(phi, xc, gamma, eta))
        Fs = np.array(Fs)

        # F(0) is the raw sharp-step IC: a central-difference gradient of a
        # true discontinuity diverges as 1/dx, a known mesh-dependent
        # discretisation artifact -- excluded from the monotonicity check
        after_ic = Fs[1:]
        diffs = np.diff(after_ic)
        max_rel_increase = float(np.max(diffs) / np.mean(after_ic)) if len(diffs) else 0.0
        is_monotone = max_rel_increase <= 1e-3
        log(f"  n={n:3d}: F(0)={Fs[0]:.4f} J/m^2 [unrelaxed step IC, discretisation "
            f"artifact]  F(final)={Fs[-1]:.4f} J/m^2  "
            f"max relative uptick after t=0: {max_rel_increase:.2e}  "
            f"monotonically non-increasing (tol 1e-3): {is_monotone}")

    # dt-halving and ramp-IC cross-checks, both against mesh10
    log("\nCross-checks (all against the mesh10 sharp-step converged profile)")
    base_dir = os.path.join(HERE, "mesh10")
    if os.path.isdir(base_dir):
        _, xc_base, phi_base, _ = converged_profile(base_dir)
        for tag, label in [("mesh10_dthalf", "half timestep"),
                            ("mesh10_ramp", "linear-ramp initial condition")]:
            case_dir = os.path.join(HERE, tag)
            if not os.path.isdir(case_dir):
                continue
            _, xc, phi, t_final = converged_profile(case_dir)
            l2, linf = l2_linf_error(phi, phi_base)
            log(f"  {tag:16s} ({label}): L2 diff vs mesh10 = {l2:.4e}, "
                f"Linf diff = {linf:.4e}")
    else:
        log("  mesh10/ not found -- run Allrun first.")

    with open(os.path.join(HERE, "convergence_table.txt"), "w") as f:
        f.write("\n".join(report_lines) + "\n")
    log("\nSaved convergence_table.txt")


if __name__ == "__main__":
    main()
