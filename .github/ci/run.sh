#!/bin/bash
# Build a minimal Cactus (cocal.th, ubuntu.cfg) with this repository as the
# COCAL arrangement, and run the COCAL_ID fixture tests. Also used by hand:
#   [JOBS=n] .github/ci/run.sh [work directory, default ./cactus-ci]
# Needs git, gcc, g++, gfortran, make, perl, python3 and OpenMPI.
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
ci=$repo/.github/ci
work=$(mkdir -p "${1:-cactus-ci}" && cd "${1:-cactus-ci}" && pwd)
cd "$work"

if [ ! -d Cactus/src ]; then
  curl -fsSL -o GetComponents https://raw.githubusercontent.com/gridaphobe/CRL/ET_2026_05/GetComponents
  chmod +x GetComponents
  ./GetComponents --parallel "$ci/cocal.th"
fi
ln -sfn "$repo" Cactus/arrangements/COCAL
{ grep -E '^[A-Za-z]+/[A-Za-z0-9_]+' "$ci/cocal.th"; echo COCAL/COCAL_ID; } > ThornList

cd Cactus
[ -d configs/cocal ] || make -j"${JOBS:-$(nproc)}" cocal-config options="$ci/ubuntu.cfg" THORNLIST="$work/ThornList" PROMPT=no
# A cached build may still list reader files that were since removed.
rm -rf configs/cocal/build/COCAL_ID
make cocal -j"${JOBS:-$(nproc)}"

# One process per fixture (test.ccl); MPI one-sided windows over shared memory.
export OMPI_MCA_osc=sm OMPI_MCA_rmaps_base_oversubscribe=1
CCTK_TESTSUITE_RUN_TESTS=COCAL_ID make cocal-testsuite PROMPT=no | tee "$work/testsuite.log"
expected=$(grep -c '^TEST ' "$repo/COCAL_ID/test/test.ccl")
passed=$(sed -n 's/.*Number of tests passed *-> *\([0-9]*\).*/\1/p' "$work/testsuite.log" | tail -1)
if [ "${passed:-0}" != "$expected" ]; then
  echo "COCAL_ID fixtures: ${passed:-0} of $expected passed" >&2
  cat TEST/cocal/COCAL_ID/*.diffs 2>/dev/null | head -100 >&2 || true
  exit 1
fi
echo "COCAL_ID fixtures: $passed of $expected passed"
