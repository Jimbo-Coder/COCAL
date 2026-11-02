# BBH reader

Set `ADMBase::initial_data = "CocalBBH"`, `COCAL_ID_PathToID`, and `Cocal`
for the initial lapse and shift. Coordinates are used as exported, without
the neutron-star radius normalization. Replace the ID paths in these examples
before use.

`COCAL_ID_bbh_fill_method = "smooth"` (default) fills both hole interiors from
their exteriors (see the [thorn README](../../README.md#black-hole-exports));
`"none"` rejects points inside the holes, for exterior-only diagnostics.
Points outside the exported patches are always an error.

## Files

For each patch `N = 1, 2, 3`:

- `bbhgrids_3D_mptN.las`: dimensions and the radial, theta and phi nodes;
- `bbhgra_3D_mptN.las`: conformal factor, lapse, shift and constants;
- `bnsgra_Kij_3D_mptN.las`: the six conformal curvature components (the BBH
  writer keeps the BNS name);
- `bin_ex_par_mptN.dat`: companion-excision and separation-node indices.

Both the 20.12 (six constants) and 23.15 (four constants) gravitational
formats are read; the curvature is 23.15. The grids are the exported nodes,
never regenerated from parameter files.

## Conventions

The holes are at `x = -dis_1` and `x = +dis_2`, `dis_1 + dis_2 = sepa`.
Patch 2 is rotated by pi about z: the x and y shift components and the xz and
yz curvature components change sign. The metric is `psi**4 delta_ij` and the
curvature the exported tensor times `psi**4`. On the hole patches psi is
fitted as r*psi from the horizon outward, and companion-excision nodes are
filled from the companion patch before fitting.

The [BBHT reader](../bbht/README.md) adds a disk to this geometry; vacuum BBH
refuses an export with `disk_state.las`.
