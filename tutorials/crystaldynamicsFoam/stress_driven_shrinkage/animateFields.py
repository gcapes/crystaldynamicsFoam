#!/usr/bin/env python3
"""Animate W_interface and n.<TRACKED_PHASE> over time into a GIF.

Run from the case directory after ./Allrun has
produced time directories. Requires numpy, matplotlib, pillow.

Usage:
    python3 animateFields.py [--phase 3] [--out W_interface_n3.gif] [--fps 2.5]
"""
import argparse
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

NY, NZ = 250, 250
Y_MM = (0.0, 200.0)
Z_MM = (0.0, 200.0)

LIST_PAT = re.compile(
    r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\((.*?)\)\s*;", re.S
)
UNIFORM_PAT = re.compile(r"internalField\s+uniform\s+([\-0-9.eE]+)\s*;")


def read_scalar_grid(path):
    text = open(path).read()
    m = LIST_PAT.search(text)
    if m:
        vals = np.array([float(x) for x in m.group(2).split()])
    else:

        um = UNIFORM_PAT.search(text)
        vals = np.full(NY * NZ, float(um.group(1)))
    return vals.reshape(NZ, NY).T  # rows=j(Y), cols=k(Z)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", type=int, default=3, help="tracked phase index (default 3)")
    ap.add_argument("--out", default=None, help="output GIF path (default W_interface_n<phase>.gif)")
    ap.add_argument("--fps", type=float, default=2.5)
    args = ap.parse_args()

    case = os.path.dirname(os.path.abspath(__file__))
    n_field = f"n.{args.phase}"
    out = args.out or os.path.join(case, f"W_interface_n{args.phase}.gif")

    times = sorted(
        (d for d in os.listdir(case) if re.match(r"^[0-9.]+$", d) and os.path.isdir(os.path.join(case, d))),
        key=float,
    )
    times = [t for t in times if os.path.exists(os.path.join(case, t, "W_interface"))]
    if not times:
        raise SystemExit(f"No written time directories with W_interface found under {case}")

    w_frames = [read_scalar_grid(os.path.join(case, t, "W_interface")) for t in times]
    n_frames = [read_scalar_grid(os.path.join(case, t, n_field)) for t in times]

    vanished_idx = next((i for i, f in enumerate(n_frames) if f.max() == 0), len(n_frames))
    calib_frames = w_frames[: max(vanished_idx, 1)]
    w_vmax = np.percentile(np.concatenate([f.ravel() for f in calib_frames]), 99.5)
    w_vmin = 0.0

    extent = [Z_MM[0], Z_MM[1], Y_MM[0], Y_MM[1]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    im1 = ax1.imshow(w_frames[0], origin="lower", cmap="inferno", vmin=w_vmin, vmax=w_vmax, extent=extent)
    ax1.set_title("W_interface")
    ax1.set_xlabel("Z (mm)")
    ax1.set_ylabel("Y (mm)")
    fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

    im2 = ax2.imshow(n_frames[0], origin="lower", cmap="viridis", vmin=0, vmax=1, extent=extent)
    ax2.set_title(n_field)
    ax2.set_xlabel("Z (mm)")
    ax2.set_ylabel("Y (mm)")
    fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

    suptitle = fig.suptitle(f"t = {times[0]}")

    def update(i):
        im1.set_data(w_frames[i])
        im2.set_data(n_frames[i])
        suptitle.set_text(f"t = {times[i]}")
        return im1, im2, suptitle

    ani = animation.FuncAnimation(fig, update, frames=len(times), interval=400, blit=False)
    fig.tight_layout()
    ani.save(out, writer=animation.PillowWriter(fps=args.fps))
    print("saved", out)


if __name__ == "__main__":
    main()
