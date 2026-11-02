# BBHT reader

Select `ADMBase::initial_data = "CocalBBHT"`, `HydroBase::initial_hydro = "Cocal"`,
and `Cocal` for the initial lapse and shift. This imports Mai Li's conformally
flat binary black hole plus circumbinary disk export. The supplied case is
unmagnetized. Coordinates and fluid quantities use `G=c=M_sun=1` directly.

## Reader organization

`src/cocal_id_bbh.f90` reads the three-patch BBH metric and, for BBHT, the
disk, which becomes a matter domain on the outer patch's grid. The metric,
hole fill and output are the ones every COCAL_ID export shares (see the
[thorn README](../../README.md)).

BHT and BBHT remain separate reader families: their metric/export formats
differ. Magnetic input belongs within each family when a supported export or
field prescription is available, using HydroBase's existing initial magnetic
field/potential selections. The current BBHT V1 path rejects Cocal magnetic
initialization because that file format contains no magnetic data. No new
magnetic prescription is assumed.

## Files and conventions

BBHT requires all [BBH metric/grid files](../bbh/README.md), plus:

- `disk_state.las`: the four unformatted sequential Fortran records described
  below. The supplied export uses little-endian data and four-byte markers.
- `diskpar.dat`: the original `bbht_disk` namelist, including `disk_rin` and
  `disk_rout`. These edges must be explicit; they are not inferred from the mesh.
- `peos_parameter.dat`: verifies that the supplied EOS has one common Gamma.
  It prints Gamma to six significant digits, so agreement is checked to that
  rounding.

The disk records are:

1. The 17-character tag `BBHT_BHT_FLUID_V1`, then three default integers
   `nr,nt,np`.
2. Double-precision `r(0:nr),theta(0:nt),phi(0:np)`, matching metric patch 3.
3. Four doubles: fitted EOS `K,Gamma`, rotation-law `logk_rot`, and Bernoulli
   constant. The latter two are not needed to import the saved fluid.
4. Double-precision `rho,Omega`, each shaped `(0:nr,0:nt,0:np)`, radial index
   fastest. `rho` is rest-mass density and `Omega=dphi/dt` about global z.

The density and Omega are the C2 spline of the saved nodes, bounded across
the disk surface as the [thorn README](../../README.md) describes, and the
density vanishes inside the exact inner edge
`disk_rin`, which lies inside a radial cell. The reader does not refill the
disk from a rotation law. Against the full-resolution export the reduced
fixture was cut from, the spline's density error in the test box is 0.01-0.1%
in the plane and up to 3.5% at z=+-1; the trilinear sampling of Mai's
`disk_sample`, used before, was 1.5-4.6% and 16%.
The saved fitted K takes precedence over the old EOS input/output K:

```
P   = K rho**Gamma
eps = K rho**(Gamma-1) / (Gamma-1)
v^i = (Omega*(-y,x,0)^i + beta^i) / alpha
W   = 1 / sqrt(1 - gamma_ij v^i v^j)
```

`eps` is specific internal energy and `v^i` is the Eulerian velocity expected
by HydroBase. Nonfinite, negative-density or non-timelike states are fatal.
Vacuum points have `rho=P=eps=v^i=0` and `W=1`.

`COCAL_ID_bbh_fill_method = "smooth"` supplies the two BH interiors by the
same automatic continuation as vacuum BBH. The disk is outside those regions.
The full Cartesian domain is evolved; this is not evolution with excision.
The fill and global patch-join limitations documented for BBH also apply here.

## Generator and tests

The collection's existing generator at
`tools/generate_parfile.py` detects `disk_state.las` and
selects BBHT. It reuses BBH grids, horizon tracking and diagnostics, and adds
Baikal, GRHayLHD, the saved EOS, rest-mass diagnostics, and a static central
refinement stack enclosing the disk. `grid.disk_levels` includes the base
level; `1` uses only the base grid for the disk. `grid.disk_radius` is the
finest disk-box half-width when refinement is enabled. The two hole
stacks retain their own resolution and optional horizon tracking.

Both `local` and `semiprod` profiles are available. These are adjustable grid
recipes, not evidence of convergence or a guarantee of fitting a given host.
The automatic static stack contains the disk while the moving hole boxes
remain inside its inner edge, including nesting and buffer clearance. Disk
boxes use the height measured from the saved disk, rather than full cubes.
Semiprod uses `dtfac=0.4`; local uses 0.25. The default bounded BBHT run is
one synchronized coarse step. Explicit checkpoint controls are supported.
BaikalVacuum is rejected for BBHT. The hybrid EOS uses the saved Gamma for
both the cold and thermal parts, reproducing a Gamma-law fluid for this case.

`bbht_mai_init.par` is an initialization-only example on a coarse smoke grid.
Replace its initial-data path before running. It includes two horizon searches,
rest mass, geometric ADM mass/angular momentum, constraints and fluid diagnostics.

Native `cocal_bbht_mai` (disk) and `cocal_bbht_mai_fill` (both filled holes)
tests use sampled original data. The fixture's source hashes and retained
indices are in `test_InitialData/bbht_mai/PROVENANCE.json`. It is for focused
reader tests and must not be used for physical evolution; sparse metric
stencils can leave uncovered exterior points. Reader failures return
`ok=.false.` with a message; the Cactus wrappers abort through `CCTK_WARN`.

Original email and attachments are preserved in `Hypatia/mai_bbh`. The source
dataset is `/data/shared/MaiLi_BBHT/unequal_light_ID` on Riemann; its preserved
copy is under `E:\phys_research\COCAL_IDs\bbht`. Collection validation records
are in `simulations/2026-10-23-cocal-bbht`.
