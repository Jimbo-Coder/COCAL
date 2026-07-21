# COCAL_ID_cpp

A pure C++ COCAL initial-data reader for Carpet-based Einstein Toolkit builds.
It supports COCAL RNS, BNS, MRNS, and BHT datasets, including piecewise
polytropic and cold tabulated EOS conversion.

The RNS and BNS implementation is a C++ port of the reduced COCAL routines in
`COCAL_ID`. The BHT raw interpolation follows PCOCAL's `coc2pri_bht_wl`. Its
automatic excision fill follows `COCAL_ID/src/Cocal/cocal_id_bht_fill.f90`:

1. Sample four exterior points along 14 radial directions.
2. Extrapolate each field to the center with fourth-order Lagrange
   interpolation in `r^2`, then average the valid directions.
3. Fill each interior ray with a cubic Hermite interpolant that matches the
   extrapolated center and the exterior value and slope at the fill radius.
4. Interpolate the tracefree extrinsic curvature separately, damp its trace
   inward, and set the excised hydro fields to vacuum.

`COCAL_ID_bht_fill_method = "smooth"` is the default and requires no auxiliary
parameters. The legacy `"IDpar"` mode requires `read_bin_AH.par` in the COCAL
ID directory and fails if that file is absent.

For MRNS data, `COCAL_ID_offset_AvecAphi = yes` evaluates the magnetic
potentials at the staggered locations used by IllinoisGRMHD. It is disabled by
default, matching `COCAL_ID`.

Run the thorn regression suite from the Cactus root with:

```sh
CCTK_TESTSUITE_RUN_TESTS=COCAL_ID_cpp \
  CCTK_TESTSUITE_RUN_PROCESSORS=1 \
  make ET_Krus-testsuite PROMPT=no
```
