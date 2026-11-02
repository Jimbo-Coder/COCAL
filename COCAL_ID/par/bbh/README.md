# BBH reader

Select `ADMBase::initial_data = "CocalBBH"`, set `COCAL_ID_PathToID`, and use
`Cocal` for the initial lapse and shift. Coordinates are used directly,
without neutron-star radius normalization. BNS, BHT, BBH and BBHT are
available in the combined development checkout.

`COCAL_ID_bbh_fill_method = "smooth"` fills both interiors automatically from
their own exterior data with the shared black-hole fill (see the
[thorn README](../../README.md#black-hole-interiors)) and no hand-set central
values. `"none"` rejects excised points for exterior-only diagnostics.
Uncovered exterior points remain fatal. The examples in this directory are
initialization inputs; replace their ID paths before use.

## Files and conventions

For each patch `N=1,2,3`, the reader requires:

- `bbhgrids_3D_mptN.las`: dimensions and actual radial/theta/phi nodes.
- `bbhgra_3D_mptN.las`: conformal factor, lapse, shift and constants.
- `bnsgra_Kij_3D_mptN.las`: six conformal curvature components.
- `bin_ex_par_mptN.dat`: companion-excision/separation-node indices.

The historical `bnsgra_Kij` name comes from the BBH writer. Both the
20.12/six-constant and 23.15/four-constant gravitational formats are
supported; curvature uses 23.15. Parameter files are not used to regenerate
the mesh, since they can describe a different grid from the saved fields.

Saved constants determine centers `x=-dis_1` and `x=+dis_2`, with
`dis_1+dis_2=sepa`. The distances need not be equal. Patch 2 is rotated by pi
about global z: x/y shift components and xz/yz curvature components change
sign, while xy does not. The physical metric is
`gamma_ij=psi**4 delta_ij`, and physical curvature is the saved tensor
multiplied by `psi**4`.

Exterior metric interpolation uses the common C2 tensor-product cubic spline
and a quintic blend over the full usable overlap. On the two hole patches,
psi is fitted as r*psi and divided by the evaluation radius. The radial fit
starts at the horizon, so the stellar centre-value replacement does not
apply. Invalid
companion-excision nodes are filled from the companion patch before fitting.
A reduced fixture with no usable overlap retains a hard switch there and
reports it. No interpolation beyond the saved finite coverage is allowed.

The spline coefficients are shared per node. The node leader reads the
exported fields and fits them; other ranks map the same storage. The BBHT disk spline is
shared the same way. The existing sharing switch can select private
per-rank storage.

## Validation and limits

The native fixtures cover equal/unequal binaries and both fill modes.
These reduced inputs are software fixtures, not physical evolutions.

Full-ID initialization, horizon recovery and bounded evolution evidence is
in the collection's `simulations/2026-10-24-bbh-bbht-validation/` record.
Spline and patch-blend constraint probes are in
`simulations/2026-10-25-patch-blend/` and
`simulations/2026-10-26-spline-hardening/`. Cartesian refinement no longer
produces the old interpolation spikes, but the exported patches still have
finite source-resolution mismatch. This is not evidence of source
convergence or long-term stability.

The [BBHT reader](../bbht/README.md) adds the supplied unmagnetized disk to
this geometry. Vacuum BBH rejects `disk_state.las` rather than silently
discarding matter. BHT remains a distinct export/fill family in the same
checkout.
