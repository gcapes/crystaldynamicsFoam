"""
Minimal OpenFOAM ASCII field / case reader used by postprocess.py.
Handles the two internalField forms crystaldynamicsFoam actually writes for a 1D
(Nx x 1 x 1) mesh: 'uniform <value>' and 'nonuniform List<scalar> N (...)'.
No external dependencies (no PyFoam / fluidfoam required).
"""
import glob
import json
import os
import re

_UNIFORM_RE = re.compile(r"internalField\s+uniform\s+([\-0-9.eE]+)\s*;")
_NONUNIFORM_RE = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(\s*(.*?)\)\s*;",
    re.S,
)


def read_scalar_field(path):
    """Return a list[float] of length Nx for a volScalarField file."""
    with open(path) as f:
        text = f.read()

    m = _NONUNIFORM_RE.search(text)
    if m:
        n = int(m.group(1))
        vals = [float(x) for x in m.group(2).split()]
        assert len(vals) == n, f"{path}: expected {n} values, got {len(vals)}"
        return vals

    m = _UNIFORM_RE.search(text)
    if m:
        # Uniform field -- caller must supply Nx to expand; return scalar.
        return float(m.group(1))

    raise ValueError(f"could not parse internalField in {path}")


def expand(value, nx):
    if isinstance(value, list):
        return value
    return [value] * nx


def list_times(case_dir):
    """Sorted list of (time_value, dirname) for numeric time directories."""
    times = []
    for name in os.listdir(case_dir):
        full = os.path.join(case_dir, name)
        if not os.path.isdir(full):
            continue
        try:
            t = float(name)
        except ValueError:
            continue
        if os.path.isfile(os.path.join(full, "n.0")):
            times.append((t, name))
    times.sort(key=lambda p: p[0])
    return times


def load_params(case_dir):
    with open(os.path.join(case_dir, "case_params.json")) as f:
        return json.load(f)


def cell_centres(params):
    eta = params["eta"]
    dx = params["dx"]
    Nx = params["Nx"]
    Lx = params["Lx"]
    return [-Lx / 2.0 + (i + 0.5) * dx for i in range(Nx)]


def find_case_dirs(root, pattern="mesh*"):
    return sorted(
        d for d in glob.glob(os.path.join(root, pattern))
        if os.path.isfile(os.path.join(d, "case_params.json"))
    )
