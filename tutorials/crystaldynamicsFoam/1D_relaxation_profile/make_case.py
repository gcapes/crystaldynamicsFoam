#!/usr/bin/env python3
"""
Generate one crystaldynamicsFoam case directory for equilibrium
double-obstacle interface profile.

Physics (from PF.H, two phases mapped to two identical materials):

    dphi/dt = M [ (4*gamma*eta/pi^2) * lap(phi) - (2*gamma/eta)*(1 - 2*phi) ]

Domain: quasi-1D strip along x, length 4*eta, 1 cell thick in y and z
(empty patches).  Interface centred at x=0.

Stability limit for the explicit phase-field source:
    dt_stable ~ dx^2 * pi^2 / (8*M*gamma*eta)

Relaxation timescale:
    tau = pi^2*eta / (4*M*gamma)

Usage:
    make_case.py OUTDIR --n-per-eta 10 [--dt-safety 0.2] [--ic step|ramp]
                          [--dt-override 1.23e-9] [--end-time-tau 12]
"""
import argparse
import json
import math
import os
import textwrap


def foam_header(cls, obj, location):
    return textwrap.dedent(f"""\
    /*--------------------------------*- C++ -*----------------------------------*\\
      =========                 |
      \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
       \\\\    /   O peration     | Website:  https://openfoam.org
        \\\\  /    A nd           | Version:  2506
         \\\\/     M anipulation  |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        format      ascii;
        class       {cls};
        location    "{location}";
        object      {obj};
    }}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    """)


def write(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(content)


def make_blockmeshdict(case_dir, Nx, dx, Lx):
    dy = dz = dx  # thin directions -- arbitrary size, cells are 'empty'
    xmin, xmax = -Lx / 2.0, Lx / 2.0
    verts = [
        (xmin, 0.0, 0.0), (xmax, 0.0, 0.0), (xmax, dy, 0.0), (xmin, dy, 0.0),
        (xmin, 0.0, dz), (xmax, 0.0, dz), (xmax, dy, dz), (xmin, dy, dz),
    ]
    vert_str = "\n".join(f"    ({v[0]:.10e} {v[1]:.10e} {v[2]:.10e}) // {i}"
                          for i, v in enumerate(verts))
    content = foam_header("dictionary", "blockMeshDict", "system") + f"""
convertToMeters 1;

vertices
(
{vert_str}
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({Nx} 1 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    xMin
    {{
        type patch;
        faces ( (0 4 7 3) );
    }}
    xMax
    {{
        type patch;
        faces ( (1 2 6 5) );
    }}
    yMin
    {{
        type empty;
        faces ( (0 1 5 4) );
    }}
    yMax
    {{
        type empty;
        faces ( (3 7 6 2) );
    }}
    zMin
    {{
        type empty;
        faces ( (0 3 2 1) );
    }}
    zMax
    {{
        type empty;
        faces ( (4 5 6 7) );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "system", "blockMeshDict"), content)


def make_controldict(case_dir, dt, end_time, write_every_steps_snapshots,
                      write_every_steps_main):
    content = foam_header("dictionary", "controlDict", "system") + f"""
libs
(
    "libROUNDSchemes.so"
);

application     crystaldynamicsFoam;

startFrom       startTime;
startTime       0;

stopAt          endTime;
endTime         {end_time:.10e};

deltaT          {dt:.10e};

// Main (AUTO_WRITE) fields are large and not needed for postprocessing here,
// so they are written only once, right at the end -- all the fields this
// case actually needs (n.0, n.1, W_interface) come from the pfSnapshots
// function object below, on its own (much finer) writeInterval.
writeControl    timeStep;
writeInterval   {write_every_steps_main};

purgeWrite      0;

writeFormat     ascii;
writePrecision  12;
writeCompression off;

timeFormat      general;
timePrecision   12;

runTimeModifiable yes;

adjustTimeStep  no;

// Not used (adjustTimeStep=no) but read unconditionally by crystaldynamicsFoam.C each step.
maxCo           1;
maxAlphaCo      1;
maxDeltaT       {dt:.10e};

functions
{{
    pfSnapshots
    {{
        type            writeObjects;
        libs            ("libutilityFunctionObjects.so");
        objects         (n.0 n.1 W_interface);
        writeControl    timeStep;
        writeInterval   {write_every_steps_snapshots};
    }}
}}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "system", "controlDict"), content)


def make_fvschemes(case_dir):
    content = foam_header("dictionary", "fvSchemes", "system") + """
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default             none;
    div(rhoPhi,U)              Gauss linear;
    "div\\(rhoPhi,n.*\\)"        Gauss ROUNDAplus01;
    div(phi,rho)                Gauss linear;
    div((Cv*rhoPhi),T)          Gauss upwind;
    div(rhoPhi,K)                Gauss upwind;
    div((sigma&U))              Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "system", "fvSchemes"), content)


def make_fvsolution(case_dir):
    content = foam_header("dictionary", "fvSolution", "system") + """
solvers
{
    "n.*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.0;
    }

    "(n.*)Final"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.0;
    }

    // Physically inert here (uniform rho, U=0 => no advection) but rhoEqn.H
    // and TEqn.H run unconditionally and need a solver entry to exist.
    rho
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.0;
    }
    rhoFinal { $rho; }

    T
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.0;
    }
    TFinal { $T; }
}

PIMPLE
{
    // nOuterCorrectors=1 so PF.H runs exactly once per timestep -- this
    // matches the single explicit-Euler-step-per-dt stability analysis in
    // the verification notes. U-solve is bypassed anyway (stressUpdate false).
    momentumPredictor   no;
    nOuterCorrectors    1;
    nCorrectors         1;
    nNonOrthogonalCorrectors 0;

    relaxU              0.5;
}

relaxationFactors
{
    equations
    {
    }
}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "system", "fvSolution"), content)


def make_grainproperties(case_dir, eta, gamma, mobility):
    # Both materials get identical properties: pairwise energy is looked up
    # per Material index with diagonal entries at 0, so the two order
    # parameters must map to two *different* materials with the *same*
    # numbers.
    material_block = """
    Cv      500;
    kappa   1.0;
    rho_0   1000;

    C0000    1.0E-1;
    C1111    1.0E-1;
    C2222    1.0E-1;

    C1010    0.1E-1;
    C2020    0.1E-1;
    C0011    0.1E-1;
    C1001    0.1E-1;
    C0022    0.1E-1;
    C2002    0.1E-1;
    C0110    0.1E-1;
    C1100    0.1E-1;
    C0101    0.1E-1;
    C2121    0.1E-1;
    C1122    0.1E-1;
    C2112    0.1E-1;
    C0220    0.1E-1;
    C2200    0.1E-1;
    C1221    0.1E-1;
    C2211    0.1E-1;
    C0202    0.1E-1;
    C1212    0.1E-1;

    N_slips  0;

    initialCRSS  16e5;
    tau_S        300e6;
    h_0          1000e6;
    r            1.4;
    k_exponent   0.0;
    gammadot_0   0.0;
"""
    content = foam_header("dictionary", "GrainProperties", "constant") + f"""
// Prescribed-stress PF-only mode: sigmaGlobal=0 so dfdp_mech is computed but
// identically zero -- isolates the interfacial physics.
stressUpdate    false;

interface_width {eta:.10e};

N_Phases        2;

N_Materials     2;

Phase2Material (0 1);

phi1 (0 0);
phi2 (0 0);
phi3 (0 0);

PF_drivingForceActive true;

n_constitutive_substeps 1;

// gamma = {gamma:.6g} J/m^2 on the only active pair; diagonals must still be
// present (Interface::Assign_Energies aborts on a missing pair) and are
// physically irrelevant here (each phase never pairs with itself).
Energy
(
(Material.0 Material.0) 0
(Material.0 Material.1) {gamma:.10e}
(Material.1 Material.1) 0
);

Mobility
(
(Material.0 Material.0) 0
(Material.0 Material.1) {mobility:.10e}
(Material.1 Material.1) 0
);

Material.0
{{
{material_block}
}}

Material.1
{{
{material_block}
}}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "constant", "GrainProperties"), content)


def scalar_field(case_dir, name, dims, internal, patches):
    boundary = "\n".join(
        f"    {p}\n    {{\n        type            {t};\n    }}"
        for p, t in patches.items()
    )
    content = foam_header("volScalarField", name, "0") + f"""
dimensions      {dims};

internalField   {internal};

boundaryField
{{
{boundary}
}}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "0", name), content)


def nonuniform_scalar(values):
    n = len(values)
    body = "\n".join(f"{v:.12e}" for v in values)
    return f"nonuniform List<scalar>\n{n}\n(\n{body}\n)"


def make_phase_fields(case_dir, Nx, dx, Lx, ic):
    xc = [-Lx / 2.0 + (i + 0.5) * dx for i in range(Nx)]

    if ic["type"] == "step":
        n0 = [1.0 if x < 0.0 else 0.0 for x in xc]
    elif ic["type"] == "ramp":
        w = ic["width"]  # full width of the ramp, centred at 0
        n0 = [min(1.0, max(0.0, 0.5 - x / w)) for x in xc]
    else:
        raise ValueError(ic["type"])
    n1 = [1.0 - v for v in n0]

    patches = {
        "xMin": "zeroGradient", "xMax": "zeroGradient",
        "yMin": "empty", "yMax": "empty", "zMin": "empty", "zMax": "empty",
    }
    scalar_field(case_dir, "n.0", "[0 0 0 0 0 0 0]", nonuniform_scalar(n0), patches)
    scalar_field(case_dir, "n.1", "[0 0 0 0 0 0 0]", nonuniform_scalar(n1), patches)


def make_U(case_dir):
    content = foam_header("volVectorField", "U", "0") + """
dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    xMin
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    xMax
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    yMin  { type empty; }
    yMax  { type empty; }
    zMin  { type empty; }
    zMax  { type empty; }
}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "0", "U"), content)


def make_T(case_dir):
    content = foam_header("volScalarField", "T", "0") + """
dimensions      [0 0 0 1 0 0 0];

internalField   uniform 300.0;

boundaryField
{
    xMin  { type zeroGradient; }
    xMax  { type zeroGradient; }
    yMin  { type empty; }
    yMax  { type empty; }
    zMin  { type empty; }
    zMax  { type empty; }
}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "0", "T"), content)


def make_sigmaGlobal(case_dir):
    content = foam_header("volTensorField", "sigmaGlobal", "0") + """
dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform (0 0 0 0 0 0 0 0 0);

boundaryField
{
    xMin  { type zeroGradient; }
    xMax  { type zeroGradient; }
    yMin  { type empty; }
    yMax  { type empty; }
    zMin  { type empty; }
    zMax  { type empty; }
}

// ************************************************************************* //
"""
    write(os.path.join(case_dir, "0", "sigmaGlobal"), content)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("outdir")
    p.add_argument("--n-per-eta", type=int, required=True,
                    help="cells across the interface width eta")
    p.add_argument("--eta", type=float, default=1.0e-6)
    p.add_argument("--gamma", type=float, default=1.0)
    p.add_argument("--mobility", type=float, default=1.0)
    p.add_argument("--dt-safety", type=float, default=0.2,
                    help="fraction of the explicit stability limit to use")
    p.add_argument("--dt-override", type=float, default=None,
                    help="use this dt exactly instead of dt-safety*dt_stable "
                         "(for the dt-halving check)")
    p.add_argument("--end-time-tau", type=float, default=3.0,
                    help="run length in units of the relaxation timescale tau")
    p.add_argument("--n-writes", type=int, default=None,
                    help="write ~this many evenly-spaced snapshots instead of "
                         "every timestep (coarser sampling, for cases that "
                         "aren't going into the per-step evolution figure)")
    p.add_argument("--ic", choices=["step", "ramp"], default="step")
    p.add_argument("--ramp-width-eta", type=float, default=3.0,
                    help="ramp width in units of eta (only used if --ic ramp)")
    args = p.parse_args()

    eta, gamma, M = args.eta, args.gamma, args.mobility
    n = args.n_per_eta
    dx = eta / n
    Lx = 4.0 * eta
    Nx = 4 * n

    dt_stable = dx ** 2 * math.pi ** 2 / (8.0 * M * gamma * eta)
    dt = args.dt_override if args.dt_override is not None else args.dt_safety * dt_stable

    tau = math.pi ** 2 * eta / (4.0 * M * gamma)
    end_time = args.end_time_tau * tau
    nsteps_total = end_time / dt

    if args.n_writes is None:
        # Default: write every single timestep, so the profile evolution
        # figure shows a smooth, unbroken transition into equilibrium.
        write_every_steps_snapshots = 1
    else:
        write_every_steps_snapshots = max(1, round(nsteps_total / args.n_writes))
    # Main (AUTO_WRITE) fields aren't needed here -- write once, near the end.
    write_every_steps_main = max(1, round(nsteps_total) + 1)

    case_dir = args.outdir
    os.makedirs(case_dir, exist_ok=True)

    make_blockmeshdict(case_dir, Nx, dx, Lx)
    make_controldict(case_dir, dt, end_time, write_every_steps_snapshots,
                      write_every_steps_main)
    make_fvschemes(case_dir)
    make_fvsolution(case_dir)
    make_grainproperties(case_dir, eta, gamma, M)

    ic = {"type": args.ic}
    if args.ic == "ramp":
        ic["width"] = args.ramp_width_eta * eta
    make_phase_fields(case_dir, Nx, dx, Lx, ic)
    make_U(case_dir)
    make_T(case_dir)
    make_sigmaGlobal(case_dir)

    open(os.path.join(case_dir, "case.foam"), "w").close()

    params = dict(
        eta=eta, gamma=gamma, mobility=M, n_per_eta=n, dx=dx, Nx=Nx, Lx=Lx,
        dt=dt, dt_stable=dt_stable, dt_safety=args.dt_safety,
        tau=tau, end_time=end_time,
        write_every_steps=write_every_steps_snapshots,
        ic=ic,
    )
    with open(os.path.join(case_dir, "case_params.json"), "w") as f:
        json.dump(params, f, indent=2)

    print(f"{case_dir}: Nx={Nx} dx={dx:.4e} dt={dt:.4e} (dt_stable={dt_stable:.4e}) "
          f"end_time={end_time:.4e} (~{args.end_time_tau:.1f} tau) "
          f"~{nsteps_total:.0f} steps, snapshot_every={write_every_steps_snapshots} "
          f"ic={args.ic}")


if __name__ == "__main__":
    main()
