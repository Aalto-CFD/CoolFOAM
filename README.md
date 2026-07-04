![OpenFOAM v14](https://img.shields.io/badge/OpenFOAM-14-brightgreen)
# CoolFOAM
The [CoolProp](https://coolprop.org/) wrapper for [OpenFOAM](https://github.com/OpenFOAM/OpenFOAM-dev).

## Overview
A fluid properties model that evaluates the thermophysical liquid properties
using the [CoolProp](https://coolprop.org) library.
It is a drop-in `liquidProperties` model, so any solver built on the standard
thermodynamics packages can use any CoolProp fluid (optionally REFPROP-,
PCSAFT-backed) without code changes.

Properties are obtained through the low-level C interface of the CoolProp
shared library, which the build downloads pre-built, so CoolProp does not need
to be compiled locally.

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

// Any CoolProp fluid name with optional explicit backend prefix
mixture
{
    // Option 1:
    H2O;

    // Option 2:
    BICUBIC&HEOS::PROPANE;

    // Option 3:
    REFPROP::IOCTANE; // requires COOLPROP_REFPROP_ROOT to be set

    // Option 4:
    PCSAFT::PROPANE;

    // See more here: https://coolprop.org/fluid_properties/PurePseudoPure.html
}
```

### Property evaluation
- Liquid-phase properties come from a pressure-temperature flash with an
  imposed liquid phase, which extends smoothly into the metastable region
  beyond saturation in the same spirit as the extrapolated built-in liquid
  correlations. If the flash fails, the saturated-liquid state at the given
  temperature is used instead.
- Vapour properties (`Cpg`, `mug`, `kappag`) use the corresponding gas-phase
  flash, and saturation properties (`pv`, `hl`, `sigma`) a quality-temperature
  flash.
- The absolute enthalpy `ha` is relative to the CoolProp reference state of the
  fluid and `hf` is `ha` at standard conditions, so the sensible energies are
  reference-independent; the absolute-energy types are therefore not meaningful
  with this model.
- The vapour diffusivity `D` is not provided by CoolProp and falls back to the
  generic API vapour mass diffusivity function with the coefficients used by
  the built-in liquids.

### Limitations
These are shared with the built-in `liquidProperties` contract:
- The equation of state is reported as incompressible (`psi = 0`, `CpMCv = 0`)
  even though the density from the pressure-temperature flash depends on
  pressure, so the compressibility terms do not see that pressure dependence.
  This is accurate away from the critical region and degrades as the critical
  point is approached.
- The saturated-liquid fallback evaluates derivative outputs (`alphav`, `Cp`)
  at the Q = 0 endpoint; a backend that refuses first derivatives exactly on
  the saturation line raises a fatal error there.
- In the region where the single-phase flash fails, every evaluation pays a
  failed pressure-temperature flash before the saturated fallback. This is the
  cost of the metastable extension and matters only for cases sitting deep
  beyond saturation.

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
