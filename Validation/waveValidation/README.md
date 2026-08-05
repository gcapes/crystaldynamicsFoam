# crystaldynamicsFoam elastic wave-speed validation

Square-pulse plane-wave benchmark on a 5 um pseudo-1-D column (2400 cells,
propagation along -y from the driven `lowerWall` patch at y = 5 um).
`Allrun` executes the case twice from the same base configuration:

| run                | pulse on `lowerWall`      | wave         | analytic speed          |
|--------------------|---------------------------|--------------|-------------------------|
| `runs/shearWave`   | (0 0 -10) m/s, tangential | shear        | c_s = sqrt(C44/rho)     |
| `runs/pressureWave`| (0 -10 0) m/s, normal     | longitudinal | c_p = sqrt(C11/rho)     |

With the distributed isotropic material (mu = 80 GPa, nu = 0.3, rho = 7800):
c_s = 3202.56 m/s and c_p = 5991.45 m/s.

`plotWaveSpeed.py` reads the probe traces (9 probes along y), detects the
wave-front first arrival at each probe (10% threshold, linearly
interpolated), fits arrival time vs position, and compares the fitted speed
with the analytic value read from `constant/GrainProperties`.  `Allrun`
returns exit status 0 only if both speeds agree within 2% (configurable via
`--tol`), so it can serve directly as an automated regression test.

Requirements: OpenFOAM (openfoam.org line, tested layout OF-10), a compiled
`crystaldynamicsFoam` in the PATH, python3 with numpy + matplotlib.

Notes
- `setFields` sets `n.0 = 1` (the initial field file is 0) -- do not skip it.
- Plasticity is disabled (`gammadot_0 0;`); the pulse stress (~250-470 MPa)
  would otherwise exceed the 16 MPa CRSS and damp the wave.
- `fvSchemes` ships with bounded vanLeer schemes so it runs on a stock
  install; swap in the ROUND library schemes for sharper history-field
  transport if compiled.
- The shear pulse is polarised along z, NOT x: `leftWall`/`rightWall` are
  `empty`, and OpenFOAM excludes empty-direction components from segregated
  vector solves (`fvMatrixSolve.C`: `validComponents[cmpt] == -1`), so an
  x-polarised pulse cannot propagate on this mesh.  `frontAndBack` uses
  zeroGradient for U (slip would clamp u_z on the z-faces and corrupt the
  boundary-face kinematics).  `plotWaveSpeed.py` cross-checks the detected
  polarisation and fails loudly if the shear run carries a longitudinal
  signal or vice versa.
- `example_output.png` shows the expected figure (generated from synthetic
  probe data at the exact analytic speeds).
