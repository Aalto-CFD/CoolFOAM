![OpenFOAM v14](https://img.shields.io/badge/OpenFOAM-14-brightgreen)
# CoolFOAM
The [CoolProp](https://coolprop.org/) wrapper for [OpenFOAM](https://github.com/OpenFOAM/OpenFOAM-dev).

## What it does
- Evaluates the thermophysical properties of a liquid or a gas through the
  [CoolProp](https://coolprop.org) library, as a drop-in `liquidProperties`
  model: any solver built on the standard thermodynamics packages can use any
  CoolProp fluid without code changes.
- Supports the optional CoolProp backends (`REFPROP::`, `PCSAFT::`,
  `BICUBIC&HEOS::`, ...) via the usual name prefix.
- Uses the low-level C interface of the CoolProp shared library, which the
  build downloads pre-built — CoolProp is never compiled locally.

## What it does not do
- Mixtures: one pure (or pseudo-pure) fluid per phase, through `pureMixture`.
- Vapour diffusivity `D` is not provided by CoolProp; it falls back to the
  generic API vapour mass diffusivity function with the built-in liquids'
  coefficients.
- Absolute-energy thermo types: `hf`/`ef` are the CoolProp reference-state
  offsets, not heats of formation, so only the sensible energies are
  meaningful.

## Usage
Load the library in `system/controlDict`:
```cpp
libs            ( "libthermophysicalProperties_aalto.so" );
```
and select the `coolprop` properties in the thermophysical dictionary:
```cpp
thermoType
{
    type        heRhoThermo;
    properties  coolprop;
    mixture     pureMixture;
    energy      sensibleInternalEnergy;
}

// Any CoolProp fluid name, with optional backend prefix:
// H2O; BICUBIC&HEOS::PROPANE; PCSAFT::PROPANE;
// REFPROP::IOCTANE (requires COOLPROP_REFPROP_ROOT to be set)
// https://coolprop.org/fluid_properties/PurePseudoPure.html
mixture
{
    H2O;

    // Optional; "liquid" (default) or "gas" -- see below
    phase           gas;
}
```

### Liquid vs gas role
The `phase` entry selects which imposed-phase flash the primary properties
(`rho`, `alphav`, `Cp`, `ha`, `ea`, `s`, `mu`, `kappa`, `psi`, `CpMCv`) come
from:

- **Liquid** (default): liquid-imposed flash, keeping the incompressible
  `liquidProperties` contract (`psi = 0`, `CpMCv = 0`, `ea = ha`).
- **Gas**: gas-imposed flash, with `psi`, `CpMCv` and `ea`/`es` reporting the
  real CoolProp compressibility, `Cp - Cv` and internal energy. Use it where
  the liquid role is not meaningful — a permanent gas such as Air at ambient
  temperature (`T >> Tc`), or the vapour phase of a two-phase VoF case.

Vapour-companion and saturation properties (`Cpg`, `mug`, `kappag`, `pv`,
`hl`, `sigma`, `pvInvert`, `Tb`, `D`) are unaffected by `phase`.

## How it works
- Primary properties come from a pressure-temperature flash with the imposed
  phase, which extends smoothly into the metastable region beyond saturation,
  in the same spirit as the extrapolated built-in liquid correlations. If the
  flash fails, the saturated state of that phase at the given temperature is
  used instead.
- Vapour-companion properties use a gas-phase flash; saturation properties a
  quality-temperature flash.
- `ha`/`ea` are relative to the CoolProp reference state of the fluid, and
  `hf`/`ef` are their values at standard conditions, so the sensible energies
  are reference-independent.

## Limitations
Shared with the built-in `liquidProperties` contract:
- **Liquid role is incompressible** (`psi = 0`, `CpMCv = 0`) even though the
  flash density does depend on pressure — accurate away from the critical
  region, degrading as it is approached. The gas role reports the real
  values and does not have this limitation.
- **Saturated fallback sits at `psat(T)`, not the requested `p`**, so beyond
  the metastable region the pressure dependence is lost and identities
  combining a primary property with the input `p` (e.g. `ha - ea` vs
  `p/rho`) hold only where the direct flash succeeds. Derivative outputs
  (`alphav`, `Cp`, gas-role `psi`/`CpMCv`) are evaluated at the saturation
  endpoint (`Q = 0` liquid, `Q = 1` gas); a backend that refuses first
  derivatives exactly on the saturation line raises a fatal error there
  (the default HEOS backend works — `Test-coolPropProperties` exercises
  this path).
- **Fallback costs a failed flash**: in the region where the single-phase
  flash fails, every evaluation pays a failed pressure-temperature flash
  first. Matters only for cases sitting deep beyond saturation.

## Compilation
CoolFOAM is normally consumed as a git submodule of Aalto's `foamSite` and
built with the rest of the site. To build it on its own, clone the repository
```sh
git clone git@github.com:Aalto-CFD/CoolFOAM.git
cd CoolFOAM
```
with OpenFOAM sourced, then compile with one of the two options below. The
build downloads the pre-built CoolProp shared library and its `CoolPropLib.h`
header into the target library directory automatically (override
`COOLPROP_VERSION` to pin a release other than the default `8.0.0`).

#### Option A: Compile into your OpenFOAM `$WM_PROJECT_SITE` (default)
```sh
wmake -all $PWD
```

#### Option B: Compile into your local user directory
```sh
FOAM_SITE_LIBBIN=$FOAM_USER_LIBBIN wmake -all $PWD
```
This places `libthermophysicalProperties_aalto.so` and the downloaded CoolProp
library under your `$FOAM_USER_LIBBIN` directory.
