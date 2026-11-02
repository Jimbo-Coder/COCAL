# COCAL_ID
A cleaned up COCAL reader for the Einstein Toolkit, Carpet based.

Parameter files for a COCAL ID can be written with the
[generator](../tools/README.md).

## How it works

Every COCAL export (RNS, MRNS, BNS) is held the same way:

- a **metric domain**: COCAL's spherical patches, one (RNS) or three (a
  binary: one per star, rotated by pi about z for the second, and the outer
  patch), with the fields psi, alpha, beta^i, K_ij/psi^4 and, for waveless
  data, h_ij;
- a **matter domain**: the fluid on its own grid, surface fitted (r/rs) as on
  COCAL's fluid grids.

At each Cactus point the metric is evaluated, the ID's matter is composed
with that same metric (the Eulerian velocity uses alpha and beta^i at the
point), and one routine writes ADMBase and HydroBase. The readers differ only
where the export formats differ.

The shift is ADMBase's contravariant beta^i. Waveless exports (RNS_WL,
MRNS_WL, BHT) store the covariant shift of the conformal metric, which the
reader raises with gt^ij; conformally flat exports need no raising.

| File | Contents |
| --- | --- |
| `src/COCAL_ID_master.F90` | Cactus entry points: read, free |
| `src/cocal_id_grid.F90` | the grid loop |
| `src/cocal_id_state.f90` | metric layout, the state at a point, EOS |
| `src/cocal_id_domain.f90` | patches, spline fit, patch blending |
| `src/cocal_id_spline.f90` | the C2 tensor-product spline |
| `src/cocal_id_eval.f90` | the state at a point, for every ID |
| `src/cocal_id_{rns,bns}.F90` | one reader per export format |
| `src/cocal_id_node_shared.*` | node-shared memory |
| `src/Cocal/` | COCAL's own routines, as exported by PCOCAL |

## Interpolation

Every field, metric and matter, is a C2 tensor-product cubic spline on its
patch: clamped in r, periodic B-splines in theta (continued over the poles)
and phi, with COCAL's fourth-order accuracy and no slope jumps at the nodes,
which the constraints would otherwise turn into errors growing like 1/h.
Grids starting at r=0 are fitted along whole diameters; the metric's centre
value is replaced by the even extrapolation of its neighbours, removing an
offset in COCAL's exports. The three patches of a binary are blended with C2
weights across their whole usable overlap, and star-patch nodes inside the
companion excision are first filled from the companion. Where an export has
no usable overlap, the reader warns and keeps a hard patch switch there.
Export files are never modified, and grids must have increasing radii and
uniform angles.

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

Matter is the spline of COCAL's Emden function on its surface-fitted grid.
As in COCAL, values at or below the vacuum threshold are vacuum.

## Node-shared storage

Node-shared arrays use the `COCAL_ID_share_node_memory` control. It defaults
to `yes`: each node's leader reads the export straight into the fields it
fits, and the spline coefficients are shared by the node's ranks. Set it to
`no` for private per-rank copies. If shared windows cannot be created at
startup the reader warns and uses private copies; a later shared allocation
failure aborts. Open MPI on Anvil needs
`OMPI_MCA_osc=sm,ucx`.

## BNS exports

Patch grids are regenerated from each patch's `rnspar_mpt*.dat`. The
coordinate map requires equal `radi` and a shared midpoint, which are checked.
Matter belongs to one star, chosen by the midplane.

### BNS velocity correction

For exports generated with the `xunit` prescription, set
`COCAL_ID::COCAL_ID_bns_xunit = yes`. This multiplies `COCAL_ID_ecc_cor_velx`
by the signed half-separation in COCAL coordinates, adding a constant x term
to each star's corotating shift; negative values give inward corrections.
The default is the radial prescription; use the prescription and coefficient
the ID was generated with.

### BNS orbital offsets

The fluid constants record in `bnsflu_3D_mpt{1,2}.las` contains
`ome, ber, radi, xcm, zcm`; missing fields are read as zero. A nonzero `zcm`
tilts the orbital rotation used for the fluid velocity.

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
