# COCAL parfile generator

[generate_parfile.py](generate_parfile.py) is one standalone Python
script shared by the main and dev branches. Generated inputs must match the
reader and evolution thorns compiled into the executable. Main provides
RNS/BNS; BBH/BBHT and per-star velocity parameters require the dev reader.
BHT and magnetic MRNS still use their dedicated parfiles.

Put the ID and run choices together in `configs` near the top, then
select an entry by its zero-based index:

```python
configs = [
    {   # 0
        "name": "bns-local",
        "id_dir": "/path/to/bns/work_area_BNS",
        "profile": "local",
        "xunit": True,
        "ecc_cor_velx": -0.11,
        "run": {"iterations": 16},
    },
    {   # 1
        "name": "baseline-semiprod",
        "id_dir": "/path/to/baseline/ID",
        "profile": "semiprod",
        "grid": {"finest_dx": 0.25, "levels": 7},
    },
]
```

```sh
python /path/to/generate_parfile.py 0
```

The selected entry chooses the profile, then overrides only the settings you
provide. The option tables below describe the available overrides.
Entries are independent; omitted values do not carry over from another run.
The shipped entry uses `/path/to/ID` as a placeholder: set that before using 0.

You can also pass a stellar ID directly to use the semiprod defaults and shared
thorn list, without selecting or inheriting anything from `configs`. A vacuum
BBH ID uses its local recipe, described below:

```sh
python /path/to/generate_parfile.py /path/to/ID
```

The optional `name` selects the output stem; otherwise it uses the ID directory
basename. Thus configuration 0 above writes `bns-local.par` and
`bns-local.json` in the current directory. Existing files are
preserved; use another directory or add `--force` to replace both outputs.
It does not build, launch or submit a simulation. The ID path, and the GRHayL
table path for tabulated IDs, are site-specific; output and checkpoint paths
are relative.

There is no separate config file or parfile template to load. The script has
three parts: editable settings at the top, parameter generation in the middle,
and a `CocalMetadata` reader class below it. It uses only the Python standard
library. Keep reusable run choices in `configs`; normal setup changes belong
there before generation, so regenerating does not lose them. Keep personal
run entries in the associated project copy and restore the generic example
before committing changes to this repository. Keep the script
with the run inputs. The JSON records its hash, selected index/name, requested
and resolved thorns, and resolved settings. Relative ID paths are relative to
the working directory; `~` is expanded. To target a directory named `0` rather
than configuration 0, pass `./0` or its absolute path.

## Select thorns

`default_thorns` near the top is the shared selection. Comment out diagnostics
you do not want by default. Use the solver selectors below to replace evolution
entries together with their settings. For a particular
run, add an `active_thorns` list inside its config; that list replaces the shared
selection for that run. Dependencies are added automatically. There are no
separate REQUIRED/OPTIONAL selection lists.

```python
    "NSSurface",          # adds SphericalSurface and interpolation
    "QuasiLocalMeasures", # adds SphericalSurface and TGRtensor
    # "Multipole",        # uncomment to add wave extraction and WeylScal4
```

For example, selecting QLM is sufficient; you need not separately list its
SphericalSurface or TGRtensor dependencies. To omit SphericalSurface entirely,
omit its consumers too: NSSurface, QLM and AHFinderDirect. Removing a diagnostic
also removes its parameter block and output requests. An explicitly selected
dependency stays active even if no other selected thorn needs it.

BNS moving grids always add VolumeIntegrals_GRMHD. RNS uses fixed boxes and can
omit both integral thorns. The CLI reports how many dependencies it added; the
JSON lists their names under `added_thorns`. `THORN_DEPENDENCIES` is the internal
dependency catalogue for this recipe, not another user selection list. It is
maintained alongside the parameter blocks when adding support for a thorn.

The executable must contain all resolved thorns. The generator neither changes
a build's thornlist nor tries to replace Cactus's own build dependencies.

## Select evolution and gauge

The shared generator owns these choices for stellar, BBH and BBHT recipes.
Saved scripts in run records are snapshots, not separate generators to maintain.
Choose hydro and spacetime independently in the same configuration:

```python
{
    "name": "bbht-grhydro-mclachlan",
    "id_dir": "/path/to/work_area_BBHT_ra0.1",
    "profile": "local",
    "hydro_thorn": "GRHydro",
    "spacetime_thorn": "ML_BSSN",
    "gauge": {
        "lapse": "1+log",
        "driver_order": 2,
        "eta": 1.0,
    },
    "run": {"iterations": 64},
}
```

| Setting | Choices | Default |
| --- | --- | --- |
| `hydro_thorn` | `GRHayLHD`, `IllinoisGRMHD`, `GRHydro` | GRHayLHD for matter; no hydro for vacuum BBH |
| `spacetime_thorn` | `Baikal`, `BaikalVacuum`, `ML_BSSN`, `ML_CCZ4`, `LeanBSSNMoL` | Baikal for matter; BaikalVacuum for vacuum BBH |

An explicit selector replaces the corresponding requested evolution thorn;
dependencies, EOS blocks, all five ADMBase evolution methods, native finite
checks, constraint variables and output requests follow the selection. Without
selectors, a custom `active_thorns` list can choose the solvers itself. Exactly
one spacetime solver and, for matter, one hydro solver must remain active.
BaikalVacuum is rejected for matter, and hydro solvers are rejected for vacuum
BBH. Existing configurations retain their defaults.

IllinoisGRMHD and GRHayLHD use the same GRHayLib EOS parameters and HydroBase
conversion cadence. GRHydro uses its own conservative variables, HydroBase
evolution method, PPM/HLLE settings and EOS_Omni's Hybrid EOS, populated with the
same cold K, Gammas, density breaks and thermal Gamma. GRHydro's tabulated-EOS
initialization/table mapping is not validated by this recipe and is rejected.
The name of a selected solver does not expand the supported COCAL input formats:
magnetized MRNS still uses its dedicated potential/staggering parfile, and the
current BBHT export has no saved magnetic field. Unmagnetized IllinoisGRMHD
cases initialize Avec/Aphi and Bvec to zero.

All five spacetime solvers can evolve vacuum; the four matter-capable choices
can be combined with each supported hydro option. These are Carpet recipes.
CarpetX evolvers use different interfaces and need their own generator recipe.
McLachlan benchmark/OpenCL variants and analytic exact-solution thorns are not
additional choices here.

`gauge` accepts:

| Key | Choices / units |
| --- | --- |
| `lapse` | `"1+log"` (default), `"harmonic"` (McLachlan only) |
| `shift` | `"gamma-driver"` |
| `driver_order` | 1 or 2 for McLachlan; fixed 2 for this Baikal build, fixed 1 for Lean |
| `advect_lapse`, `advect_shift` | Booleans, default True; fixed True for this Baikal build |
| `eta` | Nonnegative inverse evolution-length units; default `1 / source_adm_mass` |
| `lapse_floor` | Nonnegative; default 1e-15 for Baikal/McLachlan; unavailable in Lean (None) |

McLachlan maps 1+log to harmonicN=1/harmonicF=2 and harmonic slicing to
harmonicN=2/harmonicF=1. Its `evolveB` selects the driver order. For its
first-order driver, betaDriver is eta/0.75 because shiftGammaCoeff multiplies
the entire driver term; for its second-order driver betaDriver is eta.
Lean maps advection and damping to zeta_alpha/zeta_beta and eta_beta.
This Baikal build has generated gauge equations, so changing lapse family,
advection or driver order requires regenerating the thorn. Unsupported requests
fail explicitly. The isolated RK4 damping check uses the selected eta.

Slow-Start Lapse (SSL) is a separate future candidate. Etienne's
[2024 paper](https://arxiv.org/abs/2404.01137) describes a temporary damping
addition to 1+log and an open-sourced modified BaikalVacuum implementation.
The local Baikal/BaikalVacuum source has no SSL option; supporting it requires
the corresponding evolution code, not just a new generator parameter. Its
availability in a paper's code does not establish inclusion in this checkout
or validation for COCAL's initial lapse. Numerical testing is deferred.

`spacetime_options` holds additional native controls in the selected thorn's
own parameter block:

| Solver | Supported native keys |
| --- | --- |
| Baikal | `diss_strength`, `newrad_aDD_radpower` (local checkout extension) |
| BaikalVacuum | `diss_strength` |
| ML_BSSN | `epsDiss`, `conformalMethod`, `fixAdvectionTerms` |
| ML_CCZ4 | The ML_BSSN keys, plus `dampk1`, `dampk2`, `GammaShift` |
| LeanBSSNMoL | `conf_fac_floor`, `chi_gamma`, `make_aa_tracefree`, `reset_dethh` |

For example, `"spacetime_options": {"dampk1": 0.02, "dampk2": 0.0}`
sets CCZ4 damping in native units. Unspecified native controls retain the thorn
defaults; the generator does not choose an optimal CCZ4 damping prescription.
Lean uses the Dissipation thorn (fifth-order, strength 0.2); McLachlan uses
epsDiss; Baikal uses diss_strength. Equal numeric dissipation strengths do not
guarantee identical operators. Native constraint definitions and gauge
implementations also need checking when interpreting a code comparison.

The JSON records both solvers, resolved gauge and native overrides. This makes
a comparison reproducible; it does not establish numerical equivalence or
stability of a newly selected combination. Use a build containing the selected
thorns, run its `-P` parameter check, then initialize and perform a short
evolution before using a new combination for a physics study.

## ID and grid choices

Supported stellar inputs are unmagnetized BNS and RNS (CF/WL), evolved with
the selected Carpet solvers. Native `peos_parameter*.dat` or `teos_parameter*.dat`
files select the EOS type; the two BNS EOSs must agree because the hydro solver uses
one shared EOS. Tabulated IDs also need an explicit `eos_table` path to the
matching GRHayL HDF5 table in their configuration. The generator cannot infer
that table from COCAL's one-dimensional EOS export. It enables GRHayLID's
beta-equilibrium initialization and the tabulated evolution settings used by
the dedicated DD2 parfile.
Stellar radii come from the surface exports and masses from matching COCAL
diagnostics. RNS sequence selection matches the exported radius and rotation
rate. Magnetized stars still use their dedicated parfiles.

The BNS settings `compact`, `xunit` and `ecc_cor_velx` must match ID generation.
They cannot be reliably inferred from the export or folder name. The defaults
are noncompact, radial prescription, and zero correction; the command prints
the choices. For example, the new `velx110` export needs:

```python
"compact": False,
"xunit": True,
"ecc_cor_velx": -0.11,
```

A scalar `ecc_cor_velx` retains that radial/xunit prescription. For independent
Cartesian xunit coefficients, provide `ecc_cor_velx: [vx1, vx2]` and
`ecc_cor_velz: [vz1, vz2]`. These are terms in the orbital vector used before
the CO/IR/SP velocity conversion, not offsets added to the final Eulerian
velocity. Each coefficient is multiplied by `dis_cm`, retaining the existing
xunit normalization. They use signed global Cartesian components and COCAL patch 1/2 ordering, with
no additional sign reversal for patch 2. The generated runtime parameters are
`COCAL_ID_ecc_cor_velx_star[0/1]` and `COCAL_ID_ecc_cor_velz_star[0/1]`.
Both default to zero; old scalar parfiles remain valid. A scalar z value
applies the same z coefficient to both stars. Values must match the
ID generation convention; they are not inferred from the filename.

`profile` has exactly two options:

- `"local"`: eight steps, four levels, finest spacing 0.5, allocation poisoning,
  no checkpoints or 3D field output. The small domain expands with separation
  and leaves room for refinement buffers and outer-boundary stencils.
- `"semiprod"`: finest spacing 0.125, eight BNS or six RNS levels, checkpointing
  and less frequent field output. BNS uses a common box and two moving stellar
  stacks; shared depth follows the overlap of the stellar boxes. Trigger reserves
  two further central collapse levels, up to ten, activated when minimum lapse
  is below 0.1.

The `grid` dictionary has seven optional overrides; leaving it empty uses the
profile defaults:

| Setting | Meaning |
| --- | --- |
| `finest_dx` | Spacing on the finest initial level. |
| `levels` | Initial number of levels, including the base level. |
| `domain` | Three domain half-widths, aligned to coarse cells; omitted means automatic. |
| `star_margin` | Fractional surface allowance; default 0.20 (20%). |
| `star_padding` | Minimum absolute surface allowance; default 0. The larger of this and `star_margin * radius` is used. |
| `common_levels` | BNS semiprod common stack, including the base; omitted means choose from overlap. A value of 1 requests no refined central levels initially. |
| `movement_cells` | BNS displacement before moving a box, in initial finest cells; default 2. |

For example, `"grid": {"finest_dx": 0.25, "levels": 7}` gives a coarser BNS
semiprod grid with the same coarse spacing as the default. Increasing only
`finest_dx` also coarsens the outer grid and increases its timestep. To enlarge
the stellar boxes, increase `star_margin` or set `star_padding`; `common_levels`
controls the central stack independently of stellar resolution. All lengths use evolution units
(`G = c = M_sun = 1`). Local and RNS grids have no separate common stack.

`resolve_grid()` contains the whole sizing recipe:

1. Each finest box encloses its exported maximum coordinate surface radius plus
   the larger of the fractional and absolute surface allowances. BNS adds the
   box movement threshold and half a fine cell for coordinate rounding so that
   these do not consume the surface allowance.
   Half-widths round upward to a parent-cell multiple. Local boxes retain a
   minimum half-width of 12 for the small reader fixtures.
2. Moving stellar parents approximately double in size outward. Each must also
   leave conservative numerical support for its child: 12 fine cells for the
   RK4 buffers plus ghosts, and 8 parent cells for nesting distance,
   interpolation and alignment. These numbers follow this script's three ghost
   cells, four RK stages, fifth-order prolongation and `min_distance=4`.
   Carpet still builds the final union, buffers and aligned hierarchy.
3. BNS shared levels begin where the stellar boxes overlap by at least half the
   smaller box's width along the separation axis. The finest stellar level is
   kept separate. This is an efficiency choice for the default, not a Carpet
   nesting requirement. Each active shared parent encloses the next-finer
   stellar and central requests with numerical margins. Stellar requests at the
   same level can extend beyond it; Carpet uses their combined coverage.
   Wide binaries can have no refined shared levels if the overlap occurs only
   outside the available refinement hierarchy. Quasi-circular and hyperbolic
   inputs use the same construction.
4. Inactive central radii follow the stellar size scale, halving beyond the
   initial hierarchy for the collapse trigger. They do not count as initially
   active grids. The existing trigger remains centred at the origin.
5. The automatic outer domain retains the profile's minimum size and expands
   to contain the stacks with 14 coarse cells of boundary clearance. An explicit
   domain must satisfy that clearance and coarse-cell alignment. Enlarging the
   domain does not resize stellar or shared boxes. Domain size and finest
   spacing remain choices for the desired duration and accuracy.

The 20% surface allowance is a configurable baseline, not a bound on all
stellar expansion or tidal deformation. Use a larger allowance where a saved
evolution calls for it. The JSON records the initial hierarchy's stellar and
central half-widths
under `star_boxes` and `common_boxes`; only central indices below
`common_levels` are active initially (the lists start at L1). The command prints
the number of active shared refined levels, domain, finest box widths, centres,
and tracking cadence. A separate wide/close label remains only for diagnostic
sphere choices; it no longer selects a grid layout.

The `run` dictionary controls duration, timestep, output and poisoning. It is
independent of the grid overrides: a short run can retain a semiprod grid.
`walltime_hours` should match the scheduler allocation. The optional
`termination_margin_minutes` reserves time for shutdown (default 45 for stellar
semiprod, 5 for local and BBH/BBHT). Stellar `checkpoint_hours` selects periodic
walltime checkpoints (semiprod default 3, local default -1 disables periodic
checkpoints). `checkpoint_on_terminate` is independent (stellar semiprod default
True, local False). For one checkpoint at shutdown, set `checkpoint_hours=-1`
and `checkpoint_on_terminate=True`; checkpointing remains enabled and two
completed checkpoints can be retained across restarts. BBH/BBHT accept the
same explicit checkpoint settings; their default remains disabled.
For jobs limited by walltime, leave large iteration
and final-time limits; a comparison milestone need not stop the evolution.

TerminationTrigger checks walltime in ANALYSIS independently of output cadence.
Carpet performs checkpoints at synchronized meta-mode CHECKPOINT/TERMINATE
calls, so a request can wait for the next synchronization. Leave a margin for
that delay and checkpoint I/O; a scheduler hard kill cannot guarantee a final
checkpoint. Validate restart behavior with the actual executable.
Tracking runs every coarse step, independently of diagnostic/output frequency. CoM integrals,
HydroBase conversion, density-max tracking and BNS regridding use `$track_every`.
BNS matter-integral tables follow that cadence because the thorn uses one
interval for calculation and file output. Boxes move after a displacement of
`movement_cells` initial finest-grid cells: 0.25 for default semiprod, 1 locally.
RNS boxes remain at the origin.

## Evolution and diagnostics

- RK4 and a fixed Courant step: `dt = dtfac * dx` on each refinement level.
  `run["dtfac"]` defaults to 0.25 locally and 0.4 for semiprod. This is not an
  adaptive wave-speed estimate. The generator rejects the isolated Baikal
  damping violation `eta * coarse_dt >= 2.785293563`, with `eta=1/M_ADM`.
- Atmosphere/density floor `1e-12`, density ceiling `0.02`, native cold EOS,
  thermal gamma 5/3, Palenzuela1D recovery and Newman1D backup.
- AHFinderDirect searches an origin-centred horizon and exports it to the
  named SphericalSurface slot. Exterior stellar QLM spheres have radius 1.05
  times the initial maximum coordinate stellar radius and follow each CoM.
- NSSurface adds density contours in separate slots, with threshold
  `surface_density`. Shape arrays have ASCII and HDF5 output; scalar geometry
  and validity flags are saved. Keep `sf_valid`, `NSSurface::info_valid` and
  QLM validity when interpreting them.
- Global QLM and finite-volume ADM spheres must enclose the stars and leave
  four coarse cells of boundary clearance. Unsuitable spheres are omitted and
  identified in JSON. BNS semiprod retains finite ADM mass/J slots 10/11,
  12/13 and 14/15 when those spheres fit. Full-domain integrals remain, but may
  be affected by boundary errors; prefer outer QLM values for monitoring.
- Allocation/timelevel poisoning is on locally and off for semiprod. NaNChecker
  checks native GRHayL velocities at POSTSTEP and aborts on nonfinite requested
  physical fields. Carpet's separate all-storage poison scan remains off.
- Buffer zones, fifth-order spatial and second-order temporal prolongation,
  tree-based interpolation, selective storage and synchronized diagnostics are
  retained. LoopControl uses fixed tiling to avoid an initial random-search
  assertion seen locally. Semiprod enables checkpoint/recovery.

Moving BNS surfaces require the existing VolumeIntegrals_GRMHD extension
`sphericalsurface__tracks__volintegral_inside_sphere`, which is not in every
ETK checkout. Each stellar CoM moves its exterior sphere and refinement centre.
When NSSurface is enabled, the next integral reuses those same four masked
integrands to give the contour exactly the same origin, without a second mask.

## CarpetX

Set `"driver": "carpetx"` in a configuration for a CarpetX + AsterX + Z4c
parfile for COCALX's reader (Jimbo-Coder/COCALX), for RNS and BNS IDs with a
single polytrope:

```python
    {"name": "bns-carpetx", "id_dir": "/path/to/SP_K123.6_000", "driver": "carpetx"},
```

It uses the same ID metadata and grid (domain, coarse spacing, levels, boxes)
as the Carpet recipe, as CarpetX's base grid in whole blocks of 8 cells and
BoxInBox regions, and follows COCALX's examples (`COCAL_IDX/par`) for the
rest. `carpetx_thorns` near the top is its thorn selection; Cactus activates
what those thorns require. The EOS is COCAL's polytrope for the initial data and
an ideal gas of the same Gamma for the evolution, since AsterX needs a
temperature that EOSX's Hybrid EOS does not provide. `profile`, `grid`, `run`
(iterations, final_time, dtfac, out_every, out3d_every, walltime_hours,
checkpoint_hours), `gauge.eta` and the BNS velocity options apply; the Carpet
solver choices do not. Tabulated EOS, piecewise polytropes (EOSX cannot yet
set its first segment bound, 0), MRNS and black-hole IDs are refused: start
those from COCALX's examples. Its blocks are `CARPETX_BLOCKS`, one function per
thorn, like the Carpet recipe's `BLOCKS`.

## Maintenance and validation

Maintain the settings, parameter blocks and reader in this one script. The
existing committed parfiles are examples and independent inputs; changing them
no longer changes generator output. Retain each generated parfile and JSON with
its run. Updating the generator does not modify past runs or saved inputs.

`generate()` shows the overall flow for every ID: resolve the settings,
metadata, thorns, grid and run (`resolve_stellar()` or `resolve_black_holes()`),
build the shared layout (`Setup`: surface slots, integration spheres, horizon
searches), call each active thorn's block function, then render and record
provenance. The grid routine is laid out in five numbered steps, from defaults
through search radii.

Each thorn's parameters come from one function, listed in `BLOCKS` in output
order; parameter names inside it are unqualified. Stellar (RNS/BNS) and
black-hole (BBH/BBHT) IDs share each function, and where their settings differ
both are written there. To support a new thorn, add its dependencies to
`THORN_DEPENDENCIES`, write its function and add it to `BLOCKS`; only active
thorns are called, and an empty result writes no block. Shared geometry lives
in `Setup`, not in another thorn's block. `render_parfile()` handles the few
Cactus implementation aliases and formats the result.

The standalone checks reject duplicate emitted parameter names, empty values,
and invalid or qualified base names inside a block. Numerical checks cover grid
clearance, cadence, ID consistency and the isolated RK4 damping limit. They do
not parse Cactus expressions, count quotes, or scan text for `$`/`@`: those can
appear in valid expressions, environment references, escaped strings or paths.
Expression/index syntax, variable resolution, compiled thorn availability and
parameter types/ranges remain the job of `cactus_EXE -P file.par`. The generator
does not invoke an executable or claim that its own checks replace that step.

The generator preserves diagnostic slots and labels each integral separately
in the output. You can edit a generated parfile directly; regeneration requires
`--force` and replaces those edits.

The defaults preserve the established evolution choices and have short local
validation; they are not a performance optimum or a guarantee of long-run
stability. Check ID conventions, floor/ceiling activity, grid nesting and memory,
tracking, constraints/mass behaviour and restart for a new production setup.
Horizon-forming evolution and machine performance are separate validation from
startup checks. CPU layout and resource requests belong to SimFactory.

## Vacuum BBH case

A directory with `bbhgrids_3D_mpt1.las` selects the vacuum BBH recipe. It reads
all three saved metric/grid/curvature exports, companion-excision metadata and
`bbhphyseq.dat`. Saved hole centres and inner radii are used directly. A
`disk_state.las` export instead selects the separate BBHT matter recipe.
No neutron-star coordinate scaling or EOS/atmosphere settings apply to BBH.

Keep BBH choices in the same `configs` list:

```python
{
    "name": "bbh-local",
    "id_dir": "/path/to/bbhH2_test5_mr1.2_sep2.5",
    "profile": "local",
    "spacetime_thorn": "BaikalVacuum",  # or "Baikal"
    "grid": {"finest_dx": 0.0125, "levels": 5, "domain": [8, 8, 8]},
    "run": {"iterations": 64, "out_every": 16, "dtfac": 0.25},
}
```

`spacetime_thorn` selects exactly one spacetime evolution thorn. The generator switches
its parameter namespace, all five ADMBase evolution-method names, constraint
variables, NaN checks and output requests together. For regular Baikal, it
allocates TmunuBase and zeros stress-energy at the RHS stages; its scalar
reductions and NaN checks are included. BaikalVacuum omits those source terms.
Both Baikal choices describe **vacuum** here. BBHT uses a matter-capable
spacetime solver and stress-energy supplied by the fluid stack. The executable must contain the
selected thorn; it may contain both, with only one activated.

`bbh_thorns` near the top is the editable BBH selection. A per-run
`active_thorns` override retains the existing dependency/block mechanism;
`spacetime_thorn` replaces its evolution choice. The stellar selection and
stellar recipes retain their existing defaults. The executable must include
the `CocalBBH` reader to use these generated inputs.

BBH starts with two equal-size refinement stacks around the actual saved
centres. Equal stack sizes avoid a Carpet box-union/granulation failure seen
with unequal requested sizes; the physical hole radii remain independent.
Carpet constructs the union, buffer zones and proper nesting. `grid` supports
`finest_dx`, `levels`, `domain`, `track_horizons`, `movement_cells` and optional
`fine_box_radius` (one half-width for both finest boxes). The default finest half-width is the larger of 1.8
times the larger inner radius and twelve fine cells, rounded outward. A custom
width must enclose both inner radii plus four fine cells. Outer half-widths
round outward to coarse cells and must leave room around the requested stacks.

| Development choice | Grid/run settings |
| --- | --- |
| Tiny plumbing check | `finest_dx=0.05`, `levels=3`, `domain=[6,6,6]`, `iterations=4`, `out_every=4` |
| Small local smoke | `finest_dx=0.0125`, `levels=5`, `domain=[8,8,8]`, `iterations=64`, `out_every=16` |
| Finer t=0 probe | `finest_dx=0.00625`, `levels=6`, `domain=[8,8,8]`, `iterations=0`, `out_every=32` |

These numbers suit the current small-coordinate standard BBH exports, not an
arbitrary separation or mass scale. Profile defaults scale spacing with the
smaller saved inner radius: local uses R/6 and five levels; semiprod uses R/12
and six levels. Both retain the same coarse spacing and default outer
half-width 8. Both default to **four coarse steps**; semiprod chooses a finer
grid, not an unattended long evolution. A direct BBH path uses local defaults.

Carpet counts finest-level iterations: with five levels, 64 iterations are
four coarse steps. The JSON records `finest_dt`, `coarse_dt` and
`synchronized_steps`. End short runs on an output/synchronization boundary.
Boxes stay fixed by default. Set `grid.track_horizons=True` to add CarpetTracker
and follow the two AHFinderDirect surfaces with CarpetRegrid2. Horizon searches
and regrid opportunities occur every coarse synchronization, independently of
output cadence. `movement_cells` defaults to two finest cells; it controls the
displacement needed to request a move. Moving boxes snap to the coarse lattice.
The finest half-width also reserves that movement allowance beyond the four-cell
stencil margin; a custom width must provide it too. Tracking emits horizon
centroids/bounds and AMR targets, and the exterior constraint masks follow the
same AMR centres. The executable must contain CarpetTracker. Actual box positions
are quantized and lag the most recent horizon measurement.

On this Mac, moving stacks exposed a libc++ comparison issue in CarpetLib's
box-set serialization. The local seven-line comparator correction in
`CarpetLib/src/bboxset2.hh` is a separate dependency of the tracking tests;
retain that correction when reproducing them. See
`simulations/2026-10-23-cocal-bbh-tracking` in the collection for the reproducer
and run evidence. Tiny grids test the software path and cannot establish
horizon accuracy or physical convergence.

On the 16 GiB Mac, the tested tiny grid reported about 2.0 GB of active grid
storage, and the five-level local grids about 5--6 GB. The explicit six-level
unequal t=0 probe reported 18.3 GB, so use the larger-memory WSL worker for
that hierarchy. These are measured Carpet storage figures, not total-process
memory guarantees. Read the actual grid/memory report after changing spacing,
levels, box widths, domain or MPI rank count. Enlarging the domain changes
outer-grid cost without changing the finest spacing.

BBH uses AHFinderDirect's C1 Hermite geometry interpolator (`order=2`), as
recommended in that thorn's parameter documentation. Lagrange geometry
interpolation produced a missed intermediate search in the unequal smoke;
its derivatives can jump when the interpolation stencil changes. The surface
output interpolator remains Lagrange. This setting is specific to the BBH
recipe and leaves existing stellar recipes unchanged.

BBH diagnostics include two AHFinderDirect searches, QLM on the horizons and
two outer spheres, finite-volume ADM mass/J, masked exterior Hamiltonian and
momentum norms, whole-domain constraints, NaN checks, Psi4 multipoles, scalar
reductions and 1D/2D fields. There is no matter rest-mass integral for vacuum.
The exterior norm is `sqrt(integral(C^2)/integral(1))`, using coordinate volume
inside half the smallest outer half-width and excluding spheres of radius
`max(0.45, 3*R_inner)` about **both** holes. Momentum is reported by component;
whole-domain norms include the artificial fill. Masks stay fixed for fixed
boxes and follow AMR centres when tracking is enabled. NaN checks run every
step; analysis/output follows synchronized cadence.

These diagnostics being present is not numerical validation. Check
`sf_valid`, `qlm_have_valid_data` and `qlm_calc_error` before interpreting
surface quantities; invalid horizon slots can contain poison values. Outer
QLM mass on a generic sphere is not a black-hole mass. Finite-radius mass/J
from different formulas need not agree; compare matched definitions, extraction
radii and grid resolutions with the source's asymptotic values. Current smoke
tests recover horizons with the updated reader. BBH remains WIP because
constraint/source-resolution convergence and longer evolution are unvalidated.
See the collection's `simulations/2026-10-22-cocal-bbh/horizon-investigation`
record for current inputs, hashes and results.
