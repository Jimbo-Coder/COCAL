/* Node-shared storage for the COCAL_ID data.
 *
 * Every rank needs the full COCAL fluid arrays and metric spline, but they are
 * identical on all ranks. With MPI the first rank of each node allocates them
 * in an MPI-3 shared-memory window and fills them; the other ranks of the node
 * map the same memory. Without MPI, or with COCAL_ID_share_node_memory = no,
 * each rank owns a private block and is its own leader. Called from
 * cocal_id_node_shared.f90.
 *
 * Every window stays in a passive-target epoch (MPI_Win_lock_all) from
 * allocation to release, so the leader's stores are published portably by
 * MPI_Win_sync + barrier + MPI_Win_sync (MPI-3 unified memory model), not only
 * on cache-coherent x86.
 *
 * Known failure modes. An MPI that cannot create shared windows at all is
 * detected by a one-element probe window; that node then uses private copies
 * with a warning. A size limit hit later aborts with the message below (set
 * COCAL_ID_share_node_memory = no to run with private copies):
 * - Open MPI forced to osc/ucx (OMPI_MCA_osc=ucx, set by Anvil's openmpi
 *   modules): osc/ucx cannot create shared windows and returns MPI_ERR_INTERN
 *   at any size (ompi issue 11391; Anvil 4.0.6 job 21696592, also WSL 5.0.10).
 *   Allow osc/sm, e.g. OMPI_MCA_osc=sm,ucx (Anvil's RunScript does this).
 * - Shared windows live in POSIX shared memory (/dev/shm for Open MPI
 *   osc/sm), so a small /dev/shm or a container limit caps their size.
 */

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameter.h"

#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#endif

#define MAX_BLOCKS 64

static int nblocks = 0;
static void *blocks[MAX_BLOCKS];

#ifdef HAVE_CAPABILITY_MPI
static MPI_Comm node_comm = MPI_COMM_NULL;
static MPI_Win windows[MAX_BLOCKS];

static void init_node_comm(void)
{
  if (node_comm == MPI_COMM_NULL) {
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &node_comm);
    /* Report a failed window allocation ourselves, with its size. */
    MPI_Comm_set_errhandler(node_comm, MPI_ERRORS_RETURN);
  }
}

/* Collective over the node: can every rank create a shared window? Some MPI
 * setups cannot at any size (see the failure modes above); those nodes use
 * private copies instead of aborting. */
static int shared_windows_work(void)
{
  MPI_Win win;
  void *ptr;
  int ok, all, rank;
  init_node_comm();
  ok = MPI_Win_allocate_shared(sizeof(double), sizeof(double), MPI_INFO_NULL,
                               node_comm, &ptr, &win) == MPI_SUCCESS;
  MPI_Allreduce(&ok, &all, 1, MPI_INT, MPI_LAND, node_comm);
  if (all)
    MPI_Win_free(&win); /* After a partial failure a collective free could hang. */
  MPI_Comm_rank(node_comm, &rank);
  if (!all && rank == 0)
    CCTK_WARN(CCTK_WARN_ALERT,
              "COCAL_ID: this MPI cannot create node-shared windows; using "
              "private per-rank copies (more initial-data memory). See "
              "cocal_id_node_shared.c for known causes.");
  return all;
}

/* COCAL_ID::COCAL_ID_share_node_memory, and whether the MPI supports it;
 * decided once, before any rank asks whether it is the node leader. */
static int sharing(void)
{
  static int enabled = -1;
  if (enabled < 0) {
    int type;
    const CCTK_INT *value =
        CCTK_ParameterGet("COCAL_ID_share_node_memory", CCTK_THORNSTRING, &type);
    enabled = (value == NULL || *value != 0) && shared_windows_work();
  }
  return enabled;
}
#endif

/* True on the rank that reads the COCAL data and fits the spline. */
int COCAL_ID_node_leader(void)
{
#ifdef HAVE_CAPABILITY_MPI
  int rank;
  if (!sharing())
    return 1;
  init_node_comm();
  MPI_Comm_rank(node_comm, &rank);
  return rank == 0;
#else
  return 1;
#endif
}

/* Collective over the node: returns n doubles that all ranks of the node see.
 * Not zeroed: untouched pages stay non-resident until written, so large blocks
 * that are filled later do not add to the peak (cocal_id_node_shared.f90 zeroes
 * the small grid arrays). */
void *COCAL_ID_node_alloc(long n)
{
  void *ptr = NULL;

  if (nblocks == MAX_BLOCKS)
    CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID: too many node-shared blocks.");
#ifdef HAVE_CAPABILITY_MPI
  if (sharing()) {
    MPI_Aint size = COCAL_ID_node_leader() ? (MPI_Aint)n * sizeof(double) : 0;
    MPI_Aint qsize;
    int disp_unit, err;
    err = MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL,
                                  node_comm, &ptr, &windows[nblocks]);
    if (err != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int len;
      MPI_Error_string(err, msg, &len);
      CCTK_VWARN(CCTK_WARN_ABORT,
                 "COCAL_ID: cannot allocate %.2f GB of node-shared memory "
                 "(%s). Set COCAL_ID::COCAL_ID_share_node_memory = no to use "
                 "private per-rank memory instead.",
                 (double)n * sizeof(double) / 1073741824.0, msg);
    }
    MPI_Win_shared_query(windows[nblocks], 0, &qsize, &disp_unit, &ptr);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, windows[nblocks]);
  } else
#endif
  {
    ptr = malloc((size_t)n * sizeof(double));
  }
  if (ptr == NULL && n > 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "COCAL_ID: cannot allocate %.2f GB.",
               (double)n * sizeof(double) / 1073741824.0);
  blocks[nblocks++] = ptr;
  return ptr;
}

/* The leader's data must be complete before any rank of the node reads it. */
void COCAL_ID_node_barrier(void)
{
#ifdef HAVE_CAPABILITY_MPI
  int i;
  if (!sharing())
    return;
  init_node_comm();
  for (i = 0; i < nblocks; i++)
    MPI_Win_sync(windows[i]);
  MPI_Barrier(node_comm);
  for (i = 0; i < nblocks; i++)
    MPI_Win_sync(windows[i]);
#endif
}

/* Copy n doubles from the node leader to the other ranks of the node. */
void COCAL_ID_node_bcast(double *buf, int n)
{
#ifdef HAVE_CAPABILITY_MPI
  if (!sharing())
    return;
  init_node_comm();
  MPI_Bcast(buf, n, MPI_DOUBLE, 0, node_comm);
#else
  (void)buf;
  (void)n;
#endif
}

/* Collective over the node: releases every block. */
void COCAL_ID_node_free_all(void)
{
  int i;
  for (i = 0; i < nblocks; i++) {
#ifdef HAVE_CAPABILITY_MPI
    if (sharing()) {
      MPI_Win_unlock_all(windows[i]);
      MPI_Win_free(&windows[i]);
      continue;
    }
#endif
    free(blocks[i]);
  }
  nblocks = 0;
}
