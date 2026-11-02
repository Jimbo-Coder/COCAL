# BBHT reader

Set `ADMBase::initial_data = "CocalBBHT"`, `HydroBase::initial_hydro =
"Cocal"`, and `Cocal` for the initial lapse and shift. This reads Mai Li's
conformally flat binary black hole with an unmagnetized circumbinary disk, in
units G = c = M_sun = 1. The metric and hole fill are the
[BBH reader's](../bbh/README.md); the disk is a matter domain on the outer
patch's grid (`src/cocal_id_bbh.f90`). The export has no magnetic field, so
Cocal magnetic initialization is refused.

## Files

The BBH files, plus:

- `disk_state.las`: four unformatted sequential Fortran records
  (little-endian, four-byte markers):
  1. the tag `BBHT_BHT_FLUID_V1` and integers `nr, nt, np`;
  2. `r(0:nr), theta(0:nt), phi(0:np)`, the grid of metric patch 3;
  3. `K, Gamma` of the fitted polytrope, `logk_rot` and the Bernoulli
     constant (the last two unused);
  4. `rho` and `Omega = dphi/dt`, each `(0:nr,0:nt,0:np)`, r fastest;
- `diskpar.dat`: the `bbht_disk` namelist, with the edges `disk_rin` and
  `disk_rout`;
- `peos_parameter.dat`: checked for the same single Gamma (to its six printed
  digits).

## Fluid

`rho` and `Omega` are splined and bounded at the disk surface as the
[thorn README](../../README.md#interpolation) describes; `rho` is zero inside
`disk_rin`. With the saved K and Gamma:

```
P   = K rho**Gamma
eps = K rho**(Gamma-1) / (Gamma-1)
v^i = (Omega*(-y,x,0)^i + beta^i) / alpha
```

Non-finite, negative-density or non-timelike states are an error. The whole
Cartesian domain is evolved, with the holes filled as for BBH.

## Generator

`tools/generate_parfile.py` selects BBHT when the ID has `disk_state.las`. It
adds Baikal, GRHayLHD with the saved Gamma (cold and thermal), rest-mass
diagnostics and a static disk stack to the BBH grid. `grid.disk_levels`
(including the base; 1 for none) and `grid.disk_radius` (the finest disk box
half-width) set that stack; by default it holds the disk with the moving hole
boxes inside its inner edge.

`bbht_mai_init.par` is an initialization example on a coarse grid. The
`cocal_bbht_mai` and `cocal_bbht_mai_fill` tests use a reduced sample of the
export (`test_InitialData/bbht_mai/PROVENANCE.json`), not suitable for
evolution. The full export is `/data/shared/MaiLi_BBHT/unequal_light_ID` on
Riemann.
