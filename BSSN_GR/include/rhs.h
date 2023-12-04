#ifndef RHS_H
#define RHS_H

#include <time.h>

#include <cmath>
#include <iostream>

#include "block.h"
#include "derivs.h"
#include "grDef.h"
#include "mathUtils.h"
#include "parameters.h"
#include "profile_params.h"

#ifdef BSSN_ENABLE_CUDA
#include "params_cu.h"
#include "profile_gpu.h"
#include "rhs_cuda.cuh"
#endif

/**
 * @brief computes complete RHS iteratiing over all the blocks.
 *
 * This function calculates the RHS of the BSSN
 * (Baumgarte-Shapiro-Shibata-Nakamura) evolution equations based on the given
 * zipped variables, offset, bounding box information, grid size, boundary flag,
 * and optional parameters for parallel processing. Note that this will set up
 * the computations for each block and does not actually compute the RHS itself.
 *
 * @param[out] unzipVarsRHS: unzipped variables computed RHS
 * @param[in]  unzipVars: unzipped variables.
 * @param[in]  blkList: block list.
 * @param[in]  numBlocks: number of blocks.
 */
void bssnRHS(double **uzipVarsRHS, const double **uZipVars,
             const ot::Block *blkList, unsigned int numBlocks,
             unsigned int maxBlockSize);

/**
 * @brief Computes the right-hand side (RHS) of the BSSN evolution equations.
 *
 * This function calculates the RHS of the BSSN
 * (Baumgarte-Shapiro-Shibata-Nakamura) evolution equations based on the given
 * zipped variables, offset, bounding box information, grid size, boundary flag,
 * and optional parameters for parallel processing. Note that this will compute
 * the RHS for a *single* block. bssnRHS (above) will set up computation in a
 * block-wise manner and call this function
 *
 * @param[out] uzipVarsRHS Array of pointers to store the computed RHS.
 * @param[in] uZipVars Array of pointers to zipped variables.
 * @param[in] offset Offset for accessing arrays.
 * @param[in] ptmin Array representing the minimum bounding box coordinates.
 * @param[in] ptmax Array representing the maximum bounding box coordinates.
 * @param[in] sz Array representing the grid size in each dimension.
 * @param[in] bflag Boundary flag indicating the boundary conditions.
 * @param[in] threadID Identifier for the thread processing the function
 * (default is 0).
 * @param[in] maxBlockSize Maximum block size for parallel processing (default
 * is 0, indicating no specific block size, this MUST be set).
 * @param[in] blockID Identifier for the block being processed.
 */
void bssnrhs(double **uzipVarsRHS, const double **uZipVars,
             const unsigned int &offset, const double *ptmin,
             const double *ptmax, const unsigned int *sz,
             const unsigned int &bflag, const unsigned int threadID = 0,
             const unsigned int maxBlockSize = 0,
             const unsigned int blockID = 0);

void bssnrhs_sep(double **uzipVarsRHS, const double **uZipVars,
                 const unsigned int &offset, const double *ptmin,
                 const double *ptmax, const unsigned int *sz,
                 const unsigned int &bflag);

void bssn_bcs(double *f_rhs, const double *f, const double *dxf,
              const double *dyf, const double *dzf, const double *pmin,
              const double *pmax, const double f_falloff,
              const double f_asymptotic, const unsigned int *sz,
              const unsigned int &bflag);

void freeze_bcs(double *f_rhs, const unsigned int *sz,
                const unsigned int &bflag);

void fake_initial_data(double x, double y, double z, double *u);

void max_spacetime_speeds(double *const lambda1max, double *const lambda2max,
                          double *const lambda3max, const double *const alpha,
                          const double *const beta1, const double *const beta2,
                          const double *const beta3, const double *const gtd11,
                          const double *const gtd12, const double *const gtd13,
                          const double *const gtd22, const double *const gtd23,
                          const double *const gtd33, const double *const chi,
                          const unsigned int *sz);

void call_HAD_rhs();

// GRAPH AND SEPARATED PORTIONS
void bssnrhs_sep(double **unzipVarsRHS, const double **uZipVars,
                 const unsigned int &offset, const double *pmin,
                 const double *pmax, const unsigned int *sz,
                 const unsigned int &bflag);

void buildGraphVarsOnly();

void generateFullRHSGraph(int mpiRank);

void bssnrhs_sep_taskflow(double **unzipVarsRHS, const double **uZipVars,
                          const unsigned int &offset, const double *pmin,
                          const double *pmax, const unsigned int *sz,
                          const unsigned int &bflag,
                          const unsigned int threadID = 0,
                          const unsigned int maxBlockSize = 0);

void bssnRHS_TF(double **uzipVarsRHS, const double **uZipVars,
                const ot::Block *blkList, unsigned int numBlocks);

void bssnrhs_sep_taskflow_FULL(double **unzipVarsRHS, const double **uZipVars,
                               const unsigned int &offset, const double *pmin,
                               const double *pmax, const unsigned int *sz,
                               const unsigned int &bflag,
                               const unsigned int threadID = 0,
                               const unsigned int maxBlockSize = 0);

void bssnrhs_sep_taskflow_varsOnly(
    double **unzipVarsRHS, const double **uZipVars, const unsigned int &offset,
    const double *pmin, const double *pmax, const unsigned int *sz,
    const unsigned int &bflag, const unsigned int threadID = 0,
    const unsigned int maxBlockSize = 0);

void rhs_a_only(double **unzipVarsRHS, const double **uZipVars,
                const unsigned int &offset, const double *pmin,
                const double *pmax, const unsigned int *sz,
                const unsigned int &bflag, const unsigned int threadID = 0,
                const unsigned int maxBlockSize = 0);

void rhs_b_only(double **unzipVarsRHS, const double **uZipVars,
                const unsigned int &offset, const double *pmin,
                const double *pmax, const unsigned int *sz,
                const unsigned int &bflag, const unsigned int threadID = 0,
                const unsigned int maxBlockSize = 0);

void rhs_gt_only(double **unzipVarsRHS, const double **uZipVars,
                 const unsigned int &offset, const double *pmin,
                 const double *pmax, const unsigned int *sz,
                 const unsigned int &bflag, const unsigned int threadID = 0,
                 const unsigned int maxBlockSize = 0);

void rhs_chi_only(double **unzipVarsRHS, const double **uZipVars,
                  const unsigned int &offset, const double *pmin,
                  const double *pmax, const unsigned int *sz,
                  const unsigned int &bflag, const unsigned int threadID = 0,
                  const unsigned int maxBlockSize = 0);

void rhs_At_only(double **unzipVarsRHS, const double **uZipVars,
                 const unsigned int &offset, const double *pmin,
                 const double *pmax, const unsigned int *sz,
                 const unsigned int &bflag, const unsigned int threadID = 0,
                 const unsigned int maxBlockSize = 0);

void rhs_K_only(double **unzipVarsRHS, const double **uZipVars,
                const unsigned int &offset, const double *pmin,
                const double *pmax, const unsigned int *sz,
                const unsigned int &bflag, const unsigned int threadID = 0,
                const unsigned int maxBlockSize = 0);

void rhs_Gt_only(double **unzipVarsRHS, const double **uZipVars,
                 const unsigned int &offset, const double *pmin,
                 const double *pmax, const unsigned int *sz,
                 const unsigned int &bflag, const unsigned int threadID = 0,
                 const unsigned int maxBlockSize = 0);

void rhs_B_only(double **unzipVarsRHS, const double **uZipVars,
                const unsigned int &offset, const double *pmin,
                const double *pmax, const unsigned int *sz,
                const unsigned int &bflag, const unsigned int threadID = 0,
                const unsigned int maxBlockSize = 0);

#endif
