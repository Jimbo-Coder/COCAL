# COCAL_ID
A cleaned up COCAL reader for the Einstein Toolkit, Carpet based.

The [parfile generator](../tools/README.md) keeps editable run configurations
near the top of `tools/generate_parfile.py`. Select a configuration by index,
or pass an ID directory directly to use the defaults. It writes inputs without
launching a simulation.

## How it works

Every COCAL export (RNS, MRNS, BNS, BHT, BBH, BBHT) is held the same way:

- a **metric domain**: COCAL's spherical patches, one (RNS, BHT) or three (a
  binary: one per star or hole, rotated by pi about z for the second, and the
  outer patch), with the fields psi, alpha, beta^i, K_ij/psi^4 and, for
  waveless data, h_ij;
- a **matter domain**: the fluid on its own grid, surface fitted (r/rs) for
  stars as on COCAL's fluid grids, or on the gravitational grid for a torus or
  disk;
- **holes**, for black holes, filled from the exterior.

At each Cactus point the metric (or a hole's fill) is evaluated, the ID's
matter is composed with that same metric (the Eulerian velocity uses alpha
and beta^i at the point), and one routine writes ADMBase and HydroBase. The
readers differ only where the export formats differ.

The shift is ADMBase's contravariant beta^i. Waveless exports (RNS_WL,
MRNS_WL, BHT) store the covariant shift of the conformal metric, which the
reader raises with gt^ij; conformally flat exports need no raising. Readers
before November 2026 wrote the covariant shift unraised for RNS_WL and
MRNS_WL (BHT was raised), an error of about 1e-3 of the shift.

| File | Contents |
| --- | --- |
| `src/COCAL_ID_master.F90` | Cactus entry points: read, free; the same in COCALX |
| `src/cocal_id_grid.F90` | the Carpet grid loop (COCALX has its CarpetX one) |
| `src/cocal_id_state.f90` | metric layout, the state at a point, EOS |
| `src/cocal_id_domain.f90` | patches, spline fit, patch blending |
| `src/cocal_id_spline.f90` | the C2 tensor-product spline |
| `src/cocal_id_holes.f90` | black-hole interiors |
| `src/cocal_id_eval.f90` | the state at a point, for every ID |
| `src/cocal_id_{rns,bns,bht,bbh}.*` | one reader per export format |
| `src/cocal_id_node_shared.*` | node-shared memory |
| `src/Cocal/` | COCAL's own routines, as exported by PCOCAL |

## Interpolation

Every field, metric and matter, is a C2 tensor-product cubic spline on its
patch: clamped in r, periodic B-splines in theta (continued over the poles)
and phi, with COCAL's fourth-order accuracy and no slope jumps at the nodes,
which the constraints would otherwise turn into errors growing like 1/h.
Grids starting at r=0 are fitted along whole diameters; the metric's centre
value is replaced by the even extrapolation of its neighbours, removing an
offset in COCAL's exports. psi is fitted as r*psi about holes. The three
patches of a binary are blended with C2 weights across their whole usable
overlap, and star-patch nodes inside the companion excision are first filled
from the companion. Where an export has no usable overlap, the reader warns
and keeps a hard patch switch there. Export files are never modified, and
grids must have increasing radii and uniform angles.

Where a matter field's angular data are not smooth, the spline rings and
COCAL's local 4-point Lagrange interpolation does not, so there the field
blends smoothly to Lagrange in angle (the radial spline is kept). The main
case is COCAL's exported velocity-potential gradients near the poles
(`dvep_x/y`): COCAL's second-order angular differences leave them a finite
cos(2 phi) part up to the single-valued pole row. Rough nodes are found per
field from the spline coefficients: a fourth difference in theta or phi much
larger than the neighbouring second differences, and isolated, so a smooth
but under-resolved field is never marked. Smooth fields, and every field on
grids with fewer than 24 angular nodes, keep the spline. Metric fields keep
it everywhere: the constraints take its second derivatives, and their rough
spots are the binary patch seams, where Lagrange does not help.

Matter on a grid that does not fit its surface (a torus or disk on the
gravitational grid) has a kink there, across which any interpolant rings.
There the density-like field is min(spline, max(L, I*spline)), with L the
trilinear interpolant of its nodes and I that of the indicator of matter
nodes: the spline in cells of matter nodes, vacuum in cells of vacuum
nodes, at least min(spline, L) in the one layer of cells between, and
continuous everywhere. As in COCAL, values at or below the vacuum threshold
are vacuum.

## Node-shared storage

Node-shared arrays use the `COCAL_ID_share_node_memory` control. It defaults
to `yes`: each node's leader reads the export straight into the fields it
fits, and the spline coefficients are shared by the node's ranks. Set it to
`no` to use private per-rank copies. If the initial MPI capability probe
cannot create shared windows, the reader falls back to private copies with a
warning. A later large shared allocation failure aborts with guidance to
disable sharing; it does not retry privately. Open MPI on Anvil needs
`OMPI_MCA_osc=sm,ucx`.

## BNS exports

Patch grids are read from `bnsgrids_3D_mpt{1,2,3}.las` when the export
supplies them (all three or none), so the patches may have unequal radial and
angular grids; otherwise they are regenerated from each `rnspar_mpt*.dat`.
For a compactified outer export set `COCAL_ID_bns_compact = yes`: its last
radial sample, at infinity, is read but never interpolated. The coordinate
map requires equal `radi` length normalization and a shared midpoint; these
are validated constraints, distinct from the individual stellar sizes.
Matter belongs to one star, chosen by the midplane.

### BNS velocity correction

For exports generated with the `xunit` prescription, set
`COCAL_ID::COCAL_ID_bns_xunit = yes`. This multiplies `COCAL_ID_ecc_cor_velx`
by the signed half-separation in COCAL coordinates, adding a constant x term
to each star's corotating shift; negative values give inward corrections.
The default is the existing radial prescription. Use the prescription and coefficient
from the ID generation settings; neither is inferred from the directory name.

For independent coefficients use `COCAL_ID_ecc_cor_velx_star[0/1]` and
`COCAL_ID_ecc_cor_velz_star[0/1]`, in patch-1/patch-2 order. They carry their
global Cartesian signs and multiply the COCAL half-separation, `dis_cm`.
The generator accepts the pairs as `ecc_cor_velx` and `ecc_cor_velz`; scalar
x settings retain the existing radial/xunit behavior.

### BNS orbital offsets

The fluid constants record in `bnsflu_3D_mpt{1,2}.las` contains
`ome, ber, radi, xcm, zcm`. Missing `xcm` or `zcm` fields are read as zero.
Both offsets are used automatically. Nonzero `zcm` tilts the orbital rotation
used to reconstruct fluid velocities; zero `zcm` retains the existing `xcm`
calculation. The radial and `xunit` eccentricity corrections are unchanged.

## Black-hole exports

`ADMBase::initial_data = "CocalBBH"` reads conformally flat three-patch
vacuum BBH exports; see [BBH status and conventions](par/bbh/README.md).
`"CocalBBHT"` adds Mai's circumbinary disk, an unmagnetized V1 export with its
fitted single-Gamma polytrope; see [BBHT conventions](par/bbht/README.md).
`"CocalBHT"` reads the single-hole torus export, whose metric is the
one-patch waveless metric of RNS; examples are in `par/bht/`.

Each hole is filled inside its fill radius from its own patch's exterior
(`src/cocal_id_holes.f90`), in variables that keep the metric valid: log
psi^4, log alpha, the covariant shift, the tracefree curvature, h_ij and the
trace K. Centre values are extrapolated from the exterior along 14
directions, where K is taken to zero; along each ray a quintic Hermite
matches the value and first two radial derivatives of the exterior at the
fill radius (C2) and is flat at the centre. `COCAL_ID_bbh_fill_method =
"none"` rejects points inside the BBH holes instead; the BHT hole is always
filled.

See [tests and fixture provenance](test/README.md).

## Tabulated EOS utilities

`utils/make_grhayl_teos.py` converts a raw CompOSE `eos.zip`/directory to the
HDF5 tabulated-EOS schema read by GRHayL, or validates and downsamples an
existing NRPyEOS/NuLib-style HDF5 table. It requires NumPy and h5py. For
example:

```shell
uv run --with numpy --with h5py python utils/make_grhayl_teos.py \
  eos.zip DD2_GRHayL_NRPyEOS.h5 --force
```

COCAL's `teos_parameter.dat` is only a cold beta-equilibrium slice used while
mapping COCAL initial data. It is not a replacement for the hot
three-dimensional table required by GRHayLHD evolution.

`utils/validate_teos_output.py` compares COCAL_ID density and pressure output
against the COCAL TEOS slice and can write an SVG residual plot.
