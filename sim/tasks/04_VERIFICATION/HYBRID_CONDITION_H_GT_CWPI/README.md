# HYBRID_CONDITION_H_GT_CWPI

Verification case for the hybrid-code condition discussed in the Japanese
paper section on hybrid codes:

```text
h > c / omega_pi
```

The goal is not to reproduce the whole paper yet. The PDF text is not currently
machine-readable enough to extract all numerical parameters. This task isolates
the mathematical condition itself in OpenPIC.

Two subcases are provided:

```text
VALID/main.lua      h / (c / omega_pi) > 1
VIOLATED/main.lua   h / (c / omega_pi) < 1
```

Both use the same grid, field strength, and small transverse Alfvenic
perturbation. Only the background density is changed, which changes
`c / omega_pi`.

## Run

From either subcase directory:

```powershell
.\open-pic.exe .\main.lua
```

or from the build directory, pass the script path explicitly.

## Expected Use

Compare the generated grid diagnostics:

- `c_Wpi`
- `h / c_Wpi` from `h` in `lua_params.txt`
- `E`, `B`, `NP`
- growth or damping of short-wavelength artifacts

This gives a compact numerical demonstration that the hybrid condition is a
separate sanity check from the MA20/D14 square-cloud mechanism.

## Parameters

Common:

```text
nodes = 33
h     = 1 cm
B0    = 100 G
```

VALID:

```text
n0       = 1e15 cm^-3
c/Wpi    = 0.720 cm
h/c_Wpi  = 1.389
```

VIOLATED:

```text
n0       = 1e14 cm^-3
c/Wpi    = 2.277 cm
h/c_Wpi  = 0.439
```

