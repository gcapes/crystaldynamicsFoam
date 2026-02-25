import itertools

# ---------------------------------------------------------
# 1. USER INPUT: your partial stiffness entries
# ---------------------------------------------------------
# Use a dict with keys like (i,j,k,l)
C_in = {
    (0,0,0,0): 1.40e11,
    (1,1,1,1): 1.40e11,
    (2,2,2,2): 1.40e11,

    (1,0,1,0): 1.11e11,
    (2,0,2,0): 1.11e11,
    (0,0,1,1): 9.20e10,
    (1,0,0,1): 1.11e11,
    (0,0,2,2): 9.20e10,
    (2,0,0,2): 1.11e11,
    (0,1,1,0): 1.11e11,
    (1,1,0,0): 9.20e10,
    (0,1,0,1): 1.11e11,
    (2,1,2,1): 1.11e11,
    (1,1,2,2): 9.20e10,
    (2,1,1,2): 1.11e11,
    (0,2,2,0): 1.11e11,
    (2,2,0,0): 9.20e10,
    (1,2,2,1): 1.11e11,
    (2,2,1,1): 9.20e10,
    (0,2,0,2): 1.11e11,
    (1,2,1,2): 1.11e11,
}

# ---------------------------------------------------------
# 2. INITIALISE FULL 3×3×3×3 TENSOR WITH SMALL NONZERO VALUE
# ---------------------------------------------------------
# This ensures "nonzero everywhere" as you requested
C = {}
for i,j,k,l in itertools.product(range(3), repeat=4):
    C[(i,j,k,l)] = 0.0   # tiny nonzero baseline

# ---------------------------------------------------------
# 3. INSERT USER VALUES
# ---------------------------------------------------------
for key, val in C_in.items():
    C[key] = val

# ---------------------------------------------------------
# 4. ENFORCE MINOR SYMMETRY
#    C_ijkl = C_jikl = C_ijlk
# ---------------------------------------------------------
for i,j,k,l in itertools.product(range(3), repeat=4):
    v = 0.25 * (
        C[(i,j,k,l)] +
        C[(j,i,k,l)] +
        C[(i,j,l,k)] +
        C[(j,i,l,k)]
    )
    C[(i,j,k,l)] = v
    C[(j,i,k,l)] = v
    C[(i,j,l,k)] = v
    C[(j,i,l,k)] = v

# ---------------------------------------------------------
# 5. ENFORCE MAJOR SYMMETRY
#    C_ijkl = C_klij
# ---------------------------------------------------------
for i,j,k,l in itertools.product(range(3), repeat=4):
    v = 0.5 * (C[(i,j,k,l)] + C[(k,l,i,j)])
    C[(i,j,k,l)] = v
    C[(k,l,i,j)] = v

# ---------------------------------------------------------
# 6. OUTPUT IN OPENFOAM DICTIONARY FORMAT
# ---------------------------------------------------------
print("// ---- FULL 81-COMPONENT SYMMETRIC TENSOR ----")
for i,j,k,l in itertools.product(range(3), repeat=4):
    print(f"C{i}{j}{k}{l}    {C[(i,j,k,l)]:.12e};")