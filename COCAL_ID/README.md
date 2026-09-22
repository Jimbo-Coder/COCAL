# COCAL_ID
A cleaned up COCAL reader for the Einstein Toolkit, Carpet based.

## BNS velocity correction

For exports generated with the `xunit` prescription, set
`COCAL_ID::COCAL_ID_bns_xunit = yes`. This multiplies `COCAL_ID_ecc_cor_velx`
by the signed half-separation in COCAL coordinates, adding a constant x term
to each star's corotating shift; negative values give inward corrections.
The default is the existing radial prescription. Use the prescription and coefficient
from the ID generation settings; neither is inferred from the directory name.

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
