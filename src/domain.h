/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/domain.h
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "allvars.h"

#ifndef DOMAIN_H
#define DOMAIN_H

#define MASK_ACTIVE_FLAG_IN_TYPE 127
#define SET_ACTIVE_FLAG_IN_TYPE 128

extern struct local_topnode_data
{
  peanokey Size;     /*!< number of Peano-Hilbert mesh-cells represented by top-level node */
  peanokey StartKey; /*!< first Peano-Hilbert key in top-level node */
  long long Count;   /*!< counts the number of particles in this top-level node */
  double Cost;
  double SphCost;
  int Daughter; /*!< index of first daughter cell (out of 8) of top-level node */
  int Leaf;     /*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  int Parent;
  int PIndex; /*!< first particle in node */

#if defined(AMR) || defined(MODGRAV)
  int Level;
#endif
#ifdef AMR
  int MinCellLevel;
#endif
#ifdef CONSTRUCT_FOF_NGBTREE
} * topNodes, *branchNodes, *topNodes_groups; /*!< points to the root node of the top-level tree */
#else
} * topNodes, *branchNodes; /*!< points to the root node of the top-level tree */
#endif

struct domain_count_data
{
  int task;
  int count;
  int origintask;
};

extern struct domain_peano_hilbert_data
{
  peanokey key;
  int index;
#ifdef CONSTRUCT_FOF_NGBTREE
} * mp, * mp_groups;
#else
} * mp;
#endif

extern struct trans_data
{
  MyIDType ID;
  int new_task;
  int new_index;
  int wrapped;
} * trans_table;

extern int N_trans;

extern int Nbranch;

extern double fac_work, fac_load, fac_worksph;
extern double normsum_work, normsum_load, normsum_worksph;

extern double totgravcost, totpartcount, gravcost, totsphcost, sphcost;

extern struct domain_cost_data
{
  int no;
  float Work;    /*!< total "work" due to the particles stored by a leave node */
  float WorkSph; /*!< total "work" due to the particles stored by a leave node */
  int Count;     /*!< a table that gives the total number of particles held by each processor */
  int CountSph;  /*!< a table that gives the total number of SPH particles held by each processor */
#ifdef CONSTRUCT_FOF_NGBTREE
} * DomainLeaveNode, * DomainLeaveNode_groups;
#else
} * DomainLeaveNode;
#endif
/*! toGo[partner] gives the number of particles on the current task that have to go to task 'partner'
 */
extern int *toGo, *toGoSph;
extern int *toGet, *toGetSph;
extern int *list_NumPart;
extern int *list_NumGas;
extern int *list_load;
extern int *list_loadsph;
extern double *list_work;
extern double *list_worksph;

#ifdef TRACER_MC
extern int MaxNumTracer;
extern int *toGoTracer;
extern int *toGetTracer;
extern int *list_N_Tracer;
#endif
#ifdef GFM
extern int *toGoStar;
extern int *toGetStar;
extern int *list_N_star;
extern int *list_loadstar;
#endif
#ifdef BLACK_HOLES
extern int *toGoBHs;
extern int *toGetBHs;
extern int *list_NumBHs;
extern int *list_loadBHs;
#endif
#ifdef SINKS
extern int *toGoSinks;
extern int *toGetSinks;
extern int *list_NumSinks;
extern int *list_loadSinks;
#endif
#ifdef DUST_LIVE
extern int *toGoDust;
extern int *toGetDust;
extern int *list_N_dust;
extern int *list_loaddust;
#endif

void hash_table_transscribe_and_exchange(struct trans_data *trans_table);

int domain_check_for_local_refine_new(int i, MPI_Comm current_comm);
peano1D domain_double_to_int(double d);
double domain_grav_tot_costfactor(int i);
double domain_hydro_tot_costfactor(int i);
void domain_init_sum_cost(void);
void domain_printf(char *buf);
void domain_report_balance(void);
int domain_sort_load(const void *a, const void *b);
int domain_compare_count(const void *a, const void *b);
int domain_sort_task(const void *a, const void *b);
void domain_post_checks(void);
void domain_prechecks(void);
void domain_insertnode(struct local_topnode_data *treeA, struct local_topnode_data *treeB, int noA, int noB);
void domain_add_cost(struct local_topnode_data *treeA, int noA, long long count, double cost, double sphcost);
int domain_compare_count(const void *a, const void *b);
void domain_rearrange_particle_sequence(void);
void domain_combine_topleaves_to_domains(int ncpu, int ndomain);
void domain_findSplit_load_balanced(int ncpu, int ndomain);
int domain_sort_loadorigin(const void *a, const void *b);
int domain_sort_segments(const void *a, const void *b);
void domain_combine_multipledomains(void);
void domain_allocate(void);
void domain_Decomposition(void);
int domain_check_memory_bound(void);
int domain_compare_key(const void *a, const void *b);
int domain_compare_key(const void *a, const void *b);
int domain_compare_toplist(const void *a, const void *b);
double domain_particle_costfactor(int i);
int domain_countToGo(void);
int domain_decompose(void);
int domain_determineTopTree(void);

#ifdef CONSTRUCT_FOF_NGBTREE
int domain_determineTopTree_groups(void);
void domain_do_local_refine_groups(int n, int *list);
void domain_walktoptree_groups(int no);
void domain_allocate_groups(void);
void domain_free_groups(void);
void domain_allocate_lists_groups(void);
void domain_free_lists_groups(void);
#endif


void domain_exchange(void);
void domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
void domain_findExtent(void);
void domain_findSplit(int cpustart, int ncpu, int first, int last);
void domain_findSplit_balanced(int cpustart, int ncpu, int first, int last);
void domain_free(void);
void domain_shiftSplit(void);
void domain_sumCost(void);
int domain_topsplit(int node, peanokey startkey);
int domain_topsplit_local(int node, peanokey startkey, int mode);
int domain_topsplit_special(void);
int domain_compare_key(const void *a, const void *b);
int domain_check_for_local_refine(int i, MPI_Comm comm, double work);
void domain_free_trick(void);
void domain_allocate_trick(void);
int domain_recursively_combine_topTree(int start, int ncpu);
void domain_walktoptree(int no);
void domain_optimize_domain_to_task_mapping(void);
int domain_compare_count(const void *a, const void *b);
void domain_allocate_lists(void);
void domain_free_lists(void);
void domain_pack_tree_branch(int no, int parent);
int domain_unpack_tree_branch(int no, int parent);
int domain_check_for_local_refine_alt(int i, int *current_taskset);
int domain_reduce_error_flag(int flag, int *current_taskset);
void domain_do_local_refine(int n, int *list);
void domain_preserve_relevant_topnode_data(void);
void domain_find_total_cost(void);
void domain_voronoi_dynamic_update_execute(void);
void domain_prepare_voronoi_dynamic_update(void);
void domain_voronoi_dynamic_flag_particles(void);
void domain_mark_in_trans_table(int i, int task);
void domain_exchange_and_update_DC(void);
int domain_compare_connection_ID(const void *a, const void *b);
int domain_compare_local_trans_data_ID(const void *a, const void *b);
int domain_compare_recv_trans_data_ID(const void *a, const void *b);
int domain_compare_recv_trans_data_oldtask(const void *a, const void *b);

void mysort_domain(void *b, size_t n, size_t s);

#endif
