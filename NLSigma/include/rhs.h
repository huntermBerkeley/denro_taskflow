#ifndef RHS_H
#define RHS_H

#include <time.h>

#include <cmath>
#include <iostream>

#include "../../lib/taskflow/taskflow.hpp"
#include "derivs.h"
#include "mathUtils.h"
#include "nlsmUtils.h"
#include "parameters.h"

void nlsmRhs(double **uzipVarsRHS, const double **uZipVars,
             const unsigned int &offset, const double *ptmin,
             const double *ptmax, const unsigned int *sz,
             const unsigned int &bflag);

void nlsm_bcs(double *f_rhs, const double *f, const double *dxf,
              const double *dyf, const double *dzf, const double *pmin,
              const double *pmax, const double f_falloff,
              const double f_asymptotic, const unsigned int *sz,
              const unsigned int &bflag);

void fake_initial_data(double x, double y, double z, double *u);

// void nlsmRhs_blkwise();
void nlsmRHSTaskflowEntryBUILDONENTRY(double **unzipVarsRHS,
                                      const double **uZipVars,
                                      const unsigned int &offset,
                                      const double *pmin, const double *pmax,
                                      const unsigned int *sz,
                                      const unsigned int &bflag);

void nlsmRHSTaskflowFromPrebuilt(double **unzipVarsRHS, const double **uZipVars,
                                 const unsigned int &offset, const double *pmin,
                                 const double *pmax, const unsigned int *sz,
                                 const unsigned int &bflag);

void buildNLSMRHSGraph(int mpiRank);

namespace tfdendro {
struct dendrotf_rhs_data {
    double *chi;
    double *phi;
    double *phi_rhs;
    double *chi_rhs;

    unsigned int nx;
    unsigned int ny;
    unsigned int nz;

    double hx;
    double hy;
    double hz;

    int idx[3];

    unsigned int n;
    unsigned int PW;

    unsigned int BLK_SZ;

    double sigma;

    double **unzipVarsRHS;
    double **uZipVars;
    unsigned int offset;
    double *pmin;
    double *pmax;
    unsigned int *sz;
    unsigned int bflag;

    double *grad_0_chi;
    double *grad_1_chi;
    double *grad_2_chi;

    double *grad_0_phi;
    double *grad_1_phi;
    double *grad_2_phi;

    double *grad2_0_0_chi;
    double *grad2_1_1_chi;
    double *grad2_2_2_chi;
};
}  // namespace tfdendro

void make_deriv_computation_taskflow(tf::Taskflow &tf);
void make_rhs_computation_taskflow(tf::Taskflow &tf);
void make_kodiss_computation_taskflow(tf::Taskflow &tf);

#endif
