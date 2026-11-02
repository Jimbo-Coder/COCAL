# COCAL tools

[generate_parfile.py](generate_parfile.py) writes a Carpet parfile and a JSON
record for a COCAL ID. It is one script using only the Python standard
library, the same on main and dev. Main's reader supports RNS and BNS; BBH,
BBHT and the per-star velocities need dev. BHT and magnetized MRNS start from
their example parfiles.

## Usage

Runs are entries of `configs` near the top of the script:

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
python tools/generate_parfile.py 0           # configs[0]
python tools/generate_parfile.py /path/to/ID # an ID with the semiprod defaults
```

Each entry starts from its profile's defaults and overrides only what it sets.
The output is `<name>.par` and `<name>.json` in the current directory (`name`
defaults to the ID directory's name); existing files are kept unless
`--force` is given. The JSON records the script's hash, the entry, the
resolved thorns and settings. Pass `./0` for a directory named `0`.

The generator only writes inputs. Run `cactus_EXE -P file.par` to check one
against the executable.

## Thorns

`default_thorns` (stars) and `bbh_thorns` near the top are the thorn
selections; comment out a diagnostic to drop it with its parameters and
output. An entry's `active_thorns` replaces the selection for that run.
Dependencies are added from `THORN_DEPENDENCIES` (selecting QuasiLocalMeasures
adds SphericalSurface and TGRtensor). The executable must contain every
resolved thorn.

## Solvers and gauge

| Setting | Choices | Default |
| --- | --- | --- |
| `hydro_thorn` | `GRHayLHD`, `IllinoisGRMHD`, `GRHydro` | GRHayLHD; none for vacuum BBH |
| `spacetime_thorn` | `Baikal`, `BaikalVacuum`, `ML_BSSN`, `ML_CCZ4`, `LeanBSSNMoL` | Baikal; BaikalVacuum for vacuum BBH |

The solver's EOS, ADMBase evolution methods, NaN checks, constraints and
output follow the choice. BaikalVacuum is refused for matter. GRHydro uses
EOS_Omni's Hybrid EOS with the ID's polytrope; tabulated IDs need a GRHayL
solver.

`gauge`:

| Key | Choices |
| --- | --- |
| `lapse` | `"1+log"` (default), `"harmonic"` (McLachlan only) |
| `shift` | `"gamma-driver"` |
| `driver_order` | 1 or 2 for McLachlan; Baikal 2, Lean 1 |
| `advect_lapse`, `advect_shift` | default True; Baikal always True |
| `eta` | default `1/M_ADM` |
| `lapse_floor` | default 1e-15; not in Lean |

Baikal's gauge is fixed when the thorn is generated, so only `eta` and the
floor change there. `spacetime_options` passes solver parameters through:

| Solver | Keys |
| --- | --- |
| Baikal | `diss_strength`, `newrad_aDD_radpower` |
| BaikalVacuum | `diss_strength` |
| ML_BSSN | `epsDiss`, `conformalMethod`, `fixAdvectionTerms` |
| ML_CCZ4 | the ML_BSSN keys, `dampk1`, `dampk2`, `GammaShift` |
| LeanBSSNMoL | `conf_fac_floor`, `chi_gamma`, `make_aa_tracefree`, `reset_dethh` |

## Stars (RNS, BNS)

The EOS comes from the ID's `peos_parameter*.dat` or `teos_parameter*.dat`;
the two BNS stars must agree. A tabulated ID also needs `eos_table`, the
matching GRHayL HDF5 table. Radii come from the surface exports and masses
from the COCAL diagnostics.

`compact`, `xunit` and `ecc_cor_velx` must be those the BNS ID was generated
with (default noncompact, radial, 0); the command prints them. Lists
`ecc_cor_velx: [vx1, vx2]` and `ecc_cor_velz: [vz1, vz2]` give per-star
Cartesian xunit coefficients, in COCAL patch order, each multiplied by
`dis_cm`.

`profile`:

- `"local"`: 8 iterations, 4 levels, finest spacing 0.5, poisoning, no
  checkpoints or 3D output.
- `"semiprod"`: finest spacing 0.125, 8 levels (BNS) or 6 (RNS),
  checkpoints. BNS has a central stack and two stellar stacks that follow the
  stars; the collapse trigger adds up to two central levels (at most 10) when
  the lapse falls below 0.1.

`grid` overrides (lengths in units G = c = M_sun = 1):

| Key | Meaning |
| --- | --- |
| `finest_dx` | finest spacing |
| `levels` | number of levels, including the base |
| `domain` | three half-widths; default automatic |
| `star_margin` | surface allowance as a fraction of the radius, default 0.20 |
| `star_padding` | absolute surface allowance, default 0; the larger allowance is used |
| `common_levels` | BNS central stack levels, including the base; default from the stars' overlap |
| `movement_cells` | BNS displacement before a box moves, in finest cells; default 2 |

`resolve_grid()` builds the hierarchy: each finest box holds the star plus its
allowance; parents roughly double outward with room for buffers, ghosts and
prolongation; BNS levels where the stellar boxes overlap by half a box become
shared central levels; the outer domain leaves 14 coarse cells around
everything.

`run`:

| Key | Default |
| --- | --- |
| `iterations`, `final_time` | 8 local, unlimited semiprod; 1000 |
| `dtfac` | 0.25 local, 0.4 semiprod (fixed `dt = dtfac*dx`) |
| `out_every` | one coarse step; `out1d_every`, `out2d_every` multiples of it, `out3d_every` off |
| `nan_check_every` | one coarse step |
| `poison` | local only |
| `walltime_hours`, `termination_margin_minutes` | 1, 5 local; 12, 45 semiprod |
| `checkpoint_hours`, `checkpoint_on_terminate` | off local; 3 and yes semiprod |

Black-hole runs use the local margin and no checkpoints unless set. Cadences
must divide or be multiples of the coarse step. The script refuses
`eta*dt` at or above Baikal's RK4 damping limit on the coarse grid.

The stellar defaults: density floor 1e-12 and ceiling 0.02, thermal Gamma
5/3, Palenzuela1D recovery with Newman1D backup, QLM spheres at 1.05 times
each stellar radius following its centre of mass, NSSurface density contours,
and an origin AHFinderDirect search. Moving BNS spheres need VolumeIntegrals_GRMHD's
`sphericalsurface__tracks__volintegral_inside_sphere`.

## Black holes (BBH, BBHT)

`bbhgrids_3D_mpt1.las` in the ID selects the vacuum BBH recipe, and
`disk_state.las` the BBHT disk recipe. Hole centres and inner radii come from
the export.

```python
{
    "name": "bbh-local",
    "id_dir": "/path/to/bbhH2_test5_mr1.2_sep2.5",
    "profile": "local",
    "spacetime_thorn": "BaikalVacuum",
    "grid": {"finest_dx": 0.0125, "levels": 5, "domain": [8, 8, 8]},
    "run": {"iterations": 64, "out_every": 16, "dtfac": 0.25},
}
```

Both holes get equal-size refinement stacks; unequal ones trip a Carpet
box-union failure. `grid` takes `finest_dx`, `levels`, `domain`,
`movement_cells`, `fine_box_radius` (default the larger of 1.8 times the
larger inner radius and 12 fine cells) and `track_horizons` (CarpetTracker
moves the stacks with the horizons; needs CarpetTracker). Local uses spacing
R/6 and 5 levels, semiprod R/12 and 6, with R the smaller inner radius; both
default to four coarse steps.

| Check | Settings |
| --- | --- |
| plumbing | `finest_dx=0.05`, `levels=3`, `domain=[6,6,6]`, `iterations=4` |
| local smoke | `finest_dx=0.0125`, `levels=5`, `domain=[8,8,8]`, `iterations=64` |
| finer t=0 | `finest_dx=0.00625`, `levels=6`, `domain=[8,8,8]`, `iterations=0` |

Diagnostics: two AHFinderDirect searches (Hermite geometry interpolation,
`order=2`), QLM on the horizons and two outer spheres, ADM mass and J,
constraint norms outside spheres of radius `max(0.45, 3 R_inner)` about both
holes, and Psi4 multipoles. Check `sf_valid` and `qlm_have_valid_data` before
using surface values.

Moving stacks on macOS need a fix to the libc++ comparator in Carpet's
`CarpetLib/src/bboxset2.hh`.

## Tabulated EOS tables

`make_grhayl_teos.py` converts a raw CompOSE `eos.zip`/directory to the
HDF5 tabulated-EOS schema read by GRHayL, or validates and downsamples an
existing NRPyEOS/NuLib-style HDF5 table. It requires NumPy and h5py. For
example:

```shell
uv run --with numpy --with h5py python tools/make_grhayl_teos.py \
  eos.zip DD2_GRHayL_NRPyEOS.h5 --force
```

COCAL's `teos_parameter.dat` is only a cold beta-equilibrium slice used while
mapping COCAL initial data. It is not a replacement for the hot
three-dimensional table required by GRHayLHD evolution.

`validate_teos_output.py` compares COCAL_ID density and pressure output
against the COCAL TEOS slice and can write an SVG residual plot.

## Adding to the generator

`generate()` resolves the settings, ID, thorns, grid and run
(`resolve_stellar()` or `resolve_black_holes()`), builds the shared geometry
(`Setup`: surfaces, spheres, horizon searches), calls each active thorn's
function and writes the parfile and JSON. Each thorn's parameters come from
one function listed in `BLOCKS`, for stars and
black holes alike. For a new thorn, add its dependencies to
`THORN_DEPENDENCIES`, write its function and list it.
